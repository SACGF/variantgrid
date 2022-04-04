import logging
import re
from collections import defaultdict
from typing import List, Dict, Tuple, Any, Optional

import Levenshtein

from library.log_utils import log_traceback
from library.utils import is_not_none
from ontology.models import OntologyTerm, OntologyService

CodePK = Any
CodePKLookups = Dict[CodePK, str]
OntologyResults = Tuple[str, List[CodePK]]

HPO_PATTERN = re.compile(r"HP:(\d{7})$")
HPO_TYPO_PATTERN = re.compile(r"HPO:(\d{7})$")
OMIM_PATTERN = re.compile(r"OMIM:(\d+)$")
MONDO_PATTERN = re.compile(r"MONDO:(\d+)$")
HGNC_PATTERN = re.compile(r"HGNC:(\d+)$")

MIN_MATCH_LENGTH = 3
MIN_LENGTH_SINGLE_WORD_FUZZY_MATCH = 5


def load_by_id(accession, ontology_service: OntologyService) -> OntologyResults:
    pk = OntologyService.index_to_id(ontology_service, accession)
    return ontology_service, [pk]


def load_hgnc_by_id(accession) -> OntologyResults:
    return load_by_id(accession, OntologyService.HGNC)


def load_hpo_by_id(accession) -> OntologyResults:
    return load_by_id(accession, OntologyService.HPO)


def load_mondo_by_id(accession) -> OntologyResults:
    return load_by_id(accession, OntologyService.MONDO)


def load_omim_by_id(accession) -> OntologyResults:
    return load_by_id(accession, OntologyService.OMIM)


class SkipAllPhenotypeMatchException(Exception):
    pass


class PhenotypeMatcher:
    # Words which have no use matching on their own
    COMMON_WORDS = {
        'acute', 'across', 'adult', 'all', 'and', 'associated', 'auditory',
        'bad', 'bilateral', 'birth', 'blood', 'borderline', 'brain', 'brainstem',
        'can', 'carries', 'central', 'change', 'charge', 'child', 'chronic', 'close', 'comma', 'commas',
        'coned', 'cord', 'cousin', 'cousins',
        'diffused', 'deficiency', 'disorder', 'distal',
        'ear', 'exclude',
        'face', 'familial', 'floating', 'focal', 'forms', 'frequent', 'frequency', 'from', "ft 4",
        'generalized', "generalised",
        'hard', 'hearing',
        'image', 'inheritance', 'insulin',
        'joints',
        'kit',
        'large', 'lateral', 'left', 'likes', 'liver',
        'match', 'macro', 'mild', 'milena', 'moderate', 'motor', 'movements',
        'name',
        'onset', 'other',
        'parts', 'pending', 'periodic', 'person', 'pit', 'plan', 'position', 'profound', 'prolonged',
        'proximal', 'progressive',
        'range', 'raise', 'recurrent', 'right',
        'severe', 'she', 'short', 'skeletal', 'sleep', 'syndrome',
        'the', 'transient',
        'wants', 'was', 'white', 'with',
    }

    def __init__(self):
        ONTOLOGY_PK = {
            OntologyService.HPO: self._get_ontology_pks_by_term(OntologyService.HPO),
            OntologyService.MONDO: self._get_ontology_pks_by_term(OntologyService.MONDO),
            OntologyService.OMIM: self._get_omim_pks_by_term(),  # Special case
        }
        self.ontology = {}
        for ontology_service, terms_by_pk in ONTOLOGY_PK.items():
            self._break_up_terms(terms_by_pk)
            word_lookup = self._create_word_lookups(terms_by_pk)
            single_words_by_length = self._get_single_words_by_length(terms_by_pk, 5)
            self.ontology[ontology_service] = (terms_by_pk, single_words_by_length, word_lookup)

        hgnc_aliases = {}
        hgnc_names = {}
        hgnc_qs = OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC)
        for pk, name, aliases in hgnc_qs.values_list("pk", "name", "aliases"):
            for alias in aliases:
                hgnc_aliases[alias.lower()] = pk
            hgnc_names[name.lower()] = pk
        self.hgnc_records = hgnc_aliases  # Aliases first so they get overwritten
        self.hgnc_records.update(hgnc_names)  # Overwrite with assigned names

        hpo_pks = self.ontology[OntologyService.HPO][0]
        omim_pks = self.ontology[OntologyService.OMIM][0]
        special_case_lookups = self._get_special_case_lookups(hpo_pks, omim_pks, self.hgnc_records)
        self.hardcoded_lookups, self.case_insensitive_lookups, self.disease_families = special_case_lookups

    def get_matches(self, words_and_spans_subset) -> List[str]:
        words = [ws[0] for ws in words_and_spans_subset]
        text = ' '.join(words)
        lower_text = text.lower()

        ontology_term_ids: List[str] = []
        special_case_match = self._get_special_case_match(text)
        if any(special_case_match):
            ontology_term_ids = special_case_match
        else:
            if len(lower_text) < MIN_MATCH_LENGTH:
                return []

            if self._skip_word(lower_text):
                return []

            for _, (term_pks, single_words_by_length, word_lookup) in self.ontology.items():
                ontology_term_id = term_pks.get(lower_text)
                if not ontology_term_id:
                    if len(words) == 1:
                        w = words[0]
                        if len(w) >= MIN_LENGTH_SINGLE_WORD_FUZZY_MATCH:
                            ontology_term_id = self.get_id_from_single_word_fuzzy_match(single_words_by_length,
                                                                                        lower_text)
                    else:
                        lower_words = [w.lower() for w in words]
                        distance = self.calculate_match_distance(lower_words)
                        ontology_term_id = self.get_id_from_multi_word_fuzzy_match(word_lookup, lower_words, lower_text,
                                                                                   distance=distance)

                if ontology_term_id:
                    ontology_term_ids.append(ontology_term_id)

            if len(words) == 1:
                # Don't do fuzzy for genes as likely to get false positives
                if hgnc := self.hgnc_records.get(lower_text):
                    ontology_term_ids.append(hgnc)

        return ontology_term_ids

    def _get_special_case_match(self, text) -> List[str]:
        ontology_term_ids = []

        hl = self.hardcoded_lookups.get(text)
        if not hl:
            # Lookup accessions in form of: HP:0000362, OMIM:607196

            PATTERNS = [
                (HGNC_PATTERN, load_hgnc_by_id),
                (HPO_PATTERN, load_hpo_by_id),
                (HPO_TYPO_PATTERN, load_hpo_by_id),
                (MONDO_PATTERN, load_mondo_by_id),
                (OMIM_PATTERN, load_omim_by_id),
            ]

            for pattern, load_func in PATTERNS:
                if m := pattern.match(text):
                    accession = int(m.group(1))
                    hl = (load_func, accession)
                    break

        if not hl:  # Try lowercase lookups
            lowercase_text = text.lower()

            hl = self.case_insensitive_lookups.get(lowercase_text)
            if not hl:
                hl = self.disease_families.get(lowercase_text)

        if hl:
            func, arg = hl
            try:
                _, records = func(arg)
                ontology_term_ids.extend(records)

                # logging.info("Got exact: %s => %s", text, records)
            except Exception as e:
                msg = f"Error: {e}, func: {func}, arg={arg}"
                logging.error(msg)
                log_traceback()
        #            raise ValueError(msg)

        return ontology_term_ids

    @classmethod
    def _skip_word(cls, lower_text):
        """ Return true to skip a word, throws SkipAllPhenotypeMatchException to skip all.
            Only need to skip >MIN_LENGTH words as will do that later (after exact) """

        # For multi-words where you want to skip components
        SKIP_ALL = {"library prep", "to cgf", "set up", "ad pattern", "recurrent eps", "rest of"}
        if lower_text in SKIP_ALL:
            raise SkipAllPhenotypeMatchException()

        if cls._words_together(lower_text, {"TAT"}, {"non-urgent", "months", "month", "days", "week", "weeks"}):
            raise SkipAllPhenotypeMatchException()

        if cls._words_together(lower_text, {"trio"}, {"exome", "MedEx", "WES", "TS1", "father", "mother"}):
            raise SkipAllPhenotypeMatchException()

        return lower_text in cls.COMMON_WORDS

    @staticmethod
    def calculate_match_distance(words: List[str]) -> int:
        """ by default we match on 1 - however we may want to be a bit lax sometimes """
        num_ae_words = 0
        for w in words:
            num_ae_words += "ae" in w
        distance = max(1, num_ae_words)
        return distance

    @staticmethod
    def get_id_from_multi_word_fuzzy_match(lookup: CodePKLookups, words: List[str], text: str, distance: int = 1) -> Optional[CodePK]:
        potentials: CodePKLookups = dict()
        min_length = len(text) - distance
        max_length = len(text) + distance
        for w in words:
            for term, pk in lookup[w].items():
                if min_length <= len(term) <= max_length:
                    potentials[term] = pk
        if not potentials:
            return None
        return PhenotypeMatcher.get_id_from_fuzzy_match(potentials, text, distance)

    @staticmethod
    def get_id_from_single_word_fuzzy_match(single_words_by_length: Dict[int, CodePKLookups], text: str,
                                            distance: int = 1) -> Optional[CodePK]:
        text_length = len(text)
        potentials: CodePKLookups = dict()
        # Can quickly exclude words that are greater than +/- distance away
        for l in range(text_length - distance, text_length + distance + 1):
            if words := single_words_by_length.get(l):
                potentials.update(words)
        return PhenotypeMatcher.get_id_from_fuzzy_match(potentials, text, distance)

    @staticmethod
    def get_id_from_fuzzy_match(lookup: CodePKLookups, text: str, max_distance: int) -> Optional[CodePK]:
        for description, pk in lookup.items():
            distance = Levenshtein.distance(description, text)  # @UndefinedVariable
            # print("'%s' <-> '%s' distance: %d" % (description, text, distance))
            if distance <= max_distance:
                return pk
        return None

    @staticmethod
    def _words_together(text, first_words, second_words):
        first_words = {x.lower() for x in first_words}
        second_words = {x.lower() for x in second_words}

        for f in first_words:
            if f in text:
                for s in second_words:
                    if s in text:
                        return True
        return False

    @staticmethod
    def _create_word_lookups(records: CodePKLookups) -> Dict[str, Dict]:
        word_lookup = defaultdict(dict)

        for text, obj in records.items():
            for word in text.split():
                word_lookup[word][text] = obj

        return word_lookup

    @staticmethod
    def _get_single_words_by_length(records, min_length):
        """ if key is pure alpha (ie FOO not FOO,BAR or FOO1 or FOO BAR) put lookup into dict by lengths """
        words_by_length = defaultdict(dict)
        for k, v in records.items():
            if k.isalpha():
                length = len(k)
                if length >= min_length:
                    words_by_length[length][k] = v
        return words_by_length

    @staticmethod
    def _break_up_terms(hpo_pks):
        """ or could be alias """

        new_entries = {}
        for name, hpo in hpo_pks.items():
            words = name.split()
            # Make dysplasia/disease synonyms
            if len(words) >= 2:
                if words[-1] == 'dysplasia':
                    words[-1] = 'disease'
                    text = ' '.join(words)
                    if text not in hpo_pks:
                        new_entries[text] = hpo

                if words[-1] in ['dysplasia', "disease", "syndrome"]:
                    before_disease_text = ' '.join(words[:-1])
                    if before_disease_text not in hpo_pks:
                        new_entries[before_disease_text] = hpo

                # Vitamin B12 deficiency => B12 deficiency, B-12 deficiency
                if len(words) == 3 and words[0] == 'vitamin' and words[2] == 'deficiency':
                    if words[1].startswith('b'):  # Only B12 etc not "C" or "D" (need vitamin for those)
                        no_vitamin = words[1:]
                        no_vitamin_text = ' '.join(no_vitamin)
                        new_entries[no_vitamin_text] = hpo

                        b_vitamin = words[1]
                        dash_text = 'b-%s deficiency' % b_vitamin[1:]
                        new_entries[dash_text] = hpo

        hpo_pks.update(new_entries)

    @staticmethod
    def _get_ontology_pks_by_term(ontology_service: OntologyService) -> Dict[str, str]:
        ot_pks = {}
        ontology_term_qs = OntologyTerm.objects.filter(ontology_service=ontology_service)
        for pk, name, aliases in ontology_term_qs.values_list('pk', 'name', "aliases"):
            for term in filter(is_not_none, aliases + [name]):
                k = term.lower().replace(",", "")
                ot_pks[k] = pk
        return ot_pks

    @staticmethod
    def _get_omim_pks_by_term():
        """ Create entries for UNIQUE ';' separated terms
            ie "BECKWITH-WIEDEMANN SYNDROME; BWS" => "BECKWITH-WIEDEMANN" and "BECKWITH-WIEDEMANN SYNDROME" and "BWS" """

        omim_qs = OntologyTerm.objects.filter(ontology_service=OntologyService.OMIM)
        omim_pks_by_term = {}

        def break_up_dashes(omim_description, pk):
            if '-' in omim_description:
                cleaned_omim_description = omim_description.replace('-', ' ')
                omim_pks_by_term[cleaned_omim_description] = pk

        def break_up_syndromes_and_disease(omim_description, pk):
            words = omim_description.split()
            if len(words) >= 2:
                for t in ['syndrome', 'disease']:
                    try:
                        i = words.index(t)
                        if i >= 2:
                            term_before_syndrome = ' '.join(words[:i])
                            omim_pks_by_term[term_before_syndrome] = pk
                            break_up_dashes(term_before_syndrome, pk)
                    except ValueError:
                        pass

        for pk, name, aliases in omim_qs.values_list("pk", "name", "aliases"):
            for term in filter(is_not_none, aliases + [name]):
                # Remove commas, as phenotype to match will have that done also
                term = term.strip().lower().replace(",", "")
                omim_pks_by_term[term] = pk
                for split_term in term.split(";"):
                    omim_pks_by_term[split_term] = pk
                    break_up_syndromes_and_disease(split_term, pk)
                    break_up_dashes(split_term, pk)
        return omim_pks_by_term

    @staticmethod
    def _get_special_case_lookups(hpo_pks, omim_pks, gene_symbol_records) -> Tuple[Dict, Dict, Dict]:
        def load_omim_by_name(description: str) -> OntologyResults:
            omim_pk = omim_pks[description.lower()]
            return OntologyService.OMIM, [omim_pk]

        def load_hpo_by_name(hpo_name) -> OntologyResults:
            hpo_pk = hpo_pks[hpo_name.lower()]
            return OntologyService.HPO, [hpo_pk]

        def load_hpo_list_by_names(hpo_name_list) -> OntologyResults:
            hpo_list: List[CodePK] = list()
            for hpo_name in hpo_name_list:
                _, hpo = load_hpo_by_name(hpo_name)
                hpo_list.extend(hpo)
            return OntologyService.HPO, hpo_list

        def load_gene_by_name(gene_symbol: str) -> OntologyResults:
            gene_symbol_id = gene_symbol_records[gene_symbol.lower()]
            return OntologyService.HGNC, [gene_symbol_id]

        def load_genes_by_name(gene_symbols_list: List[str]) -> OntologyResults:
            genes_list = []
            for gene_symbol in gene_symbols_list:
                _, genes = load_gene_by_name(gene_symbol)
                genes_list.extend(genes)

            return OntologyService.HGNC, genes_list

        omim_qs = OntologyTerm.objects.filter(ontology_service=OntologyService.OMIM)

        def load_omim_pks_containing_name(name: str) -> OntologyResults:
            return OntologyService.OMIM, omim_qs.filter(name__icontains=name).values_list("pk", flat=True)

        def load_omim_pks_containing_alias_name(name) -> OntologyResults:
            return OntologyService.OMIM, omim_qs.filter(aliases__icontains=name).values_list("pk", flat=True)

        ABSENT_FOREARM = (load_hpo_by_name, 'absent forearm')
        ABNORMAL_BRAIN = (load_hpo_by_name, "Abnormality of brain morphology")
        ABNORMALITY_OF_LIMBS = (load_hpo_by_name, 'Abnormality of limbs')
        ARYLSULFATASE_A_DEFICIENCY = (load_omim_by_name, "ARYLSULFATASE A DEFICIENCY")
        AUTISTIC = (load_hpo_by_id, 729)
        BULLS_EYE_MACULOPATHY = (load_hpo_by_name, "bull's eye maculopathy")
        FATTY_ACID_DISORDER = (load_hpo_by_id, 4359)
        HUS = (load_hpo_by_name, "Hemolytic-uremic syndrome")
        MITO_DEFICIENCY = (load_omim_by_name, "MITOCHONDRIAL COMPLEX I DEFICIENCY")
        DEVELOPMENTAL_DELAY = (load_hpo_by_id, 1263)
        ELEVATED_CK = (load_hpo_by_id, 30234)  # Highly elevated CK
        KETOSIS = (load_hpo_by_name, "Ketosis")
        PIERRE_ROBIN = (load_hpo_by_name, "Pierre-Robin sequence")
        PAVM = (load_hpo_by_name, "Pulmonary arteriovenous malformation")
        CMS = (load_hpo_by_name, "Fatigable weakness")
        HYDROPS_FETALIS = (load_hpo_by_name, "Nonimmune hydrops fetalis")
        AFEBRILE = (load_hpo_by_id, 7359)  # Focal seizures, afebrile (HP:0040168) is obsolete, links to "Focal-onset seizure"
        FEBRILE_SEIZURES = (load_hpo_by_id, 11171)
        GEFS = FEBRILE_SEIZURES  # GEFS+ is a multi-type OMIM disease (febrile seizures links to all)
        PIG_GENES = (load_genes_by_name, ['PIG' + i for i in 'ABCFGHKLMNOPQSTUVWXYZ'])
        GLYCOGEN_STORAGE_DISEASE = (load_omim_pks_containing_name, "glycogen storage disease")
        HIGH_TSH = (load_hpo_by_id, 2925)  # Increased thyroid-stimulating hormone level
        PARKINSONISM = (load_hpo_by_name, 'Parkinsonism')
        DIBETES_TYPE_1 = (load_hpo_by_name, 'Type I diabetes mellitus')
        DYSMORPHIC_FACE = (load_hpo_by_name, "Abnormal facial shape")
        COGNITIVE_IMPAIRMENT = (load_hpo_by_name, "Cognitive impairment")  # There is also 'Specific learning disability' but this is different I think
        FACIAL_DYSMORPHISM = (load_hpo_by_id, 1999)
        HEARING_IMPAIRMENT = (load_hpo_by_name, "Hearing impairment")

        HARDCODED_LOOKUPS = {
            'aHUS': HUS,
            "ALL": (load_hpo_by_name, "Acute lymphoblastic leukemia"),
            # AML fix until we get new HPO data - see https://github.com/obophenotype/human-phenotype-ontology/issues/4236
            "AML": (load_hpo_by_name, "Acute myeloid leukemia"),
            "ADPCKD": (load_omim_by_id, 600273),  # Autosomal dominant polycystic kidney disease
            "AVSD": (load_hpo_by_name, "Atrioventricular canal defect"),  # aka Atrioventricular septal defect
            "BCC": (load_hpo_by_name, "Basal cell carcinoma"),
            "BrCa": (load_omim_by_id, 114480),  # BREAST CANCER
            "CHD": (load_hpo_by_name, "Abnormal heart morphology"),
            "CMS": CMS,
            "DD":  DEVELOPMENTAL_DELAY,
            "FAOD": FATTY_ACID_DISORDER,  # Fatty Acid Oxidation Disorders
            "FSGS": (load_hpo_by_name, "focal segmental glomerulosclerosis"),
            "FTT": (load_hpo_by_name, "Failure to thrive"),
            "GAII": (load_omim_by_name, "GLUTARIC ACIDURIA II"),
            "GEFS": GEFS,
            "GEFS+": GEFS,
            "GSD": GLYCOGEN_STORAGE_DISEASE,
            "GTOP": (load_hpo_by_name, "Spontaneous abortion"),  # Genetic Termination of Pregnancy
            "HCM": (load_hpo_by_name, "Concentric hypertrophic cardiomyopathy"),
            "HL": (load_hpo_by_name, "Hodgkin lymphoma"),
            'HUS': HUS,
            "IBD": (load_omim_by_id, 266600),  # IBD1
            "ID": (load_hpo_by_name, 'intellectual disability'),
            "LGA": (load_hpo_by_name, "Large for gestational age"),
            "LQTS": (load_hpo_by_id, 31547),  # Long QT syndrome
            "MM": (load_hpo_by_name, 'Multiple myeloma'),
            "NCS": (load_hpo_by_id, 12668),  # "Neurocardiogenic syncope" aka Vasovagal syncope
            "PCKD": (load_hpo_by_name, "Polycystic kidney dysplasia"),
            "PV": (load_omim_by_id, 263300),  # POLYCYTHEMIA VERA; PV
            "SCID": (load_hpo_by_name, "Severe combined immunodeficiency"),
            'SMA': (load_hpo_by_name, "spinal muscular atrophy"),
            "SNA12": (load_gene_by_name, "SNAI2"),  # Common misspelling
            "SUDEP": (load_hpo_list_by_names, ["Sudden death", "Epilepsy"]),
            "VSD": (load_hpo_by_name, "Ventricular septal defect"),
        }

        CASE_INSENSITIVE_LOOKUPS = {
            "aarskog": (load_omim_by_name, "AARSKOG-SCOTT SYNDROME"),
            "abdo pain": (load_hpo_by_name, "Abdominal pain"),
            "abnormal mri brain": ABNORMAL_BRAIN,
            "aching limbs": (load_hpo_by_name, "Limb pain"),
            "adenosine phosphoribosyl transferase deficiencies": (load_omim_by_id, 614723),
            "agenesis cc": (load_hpo_by_id, 1274),  # Agenesis of corpus callosum
            "afebrile seizures": AFEBRILE,
            "autistic features": AUTISTIC,
            "autistic": AUTISTIC,
            "behaviour problems": (load_hpo_by_name, "Behavioral abnormality"),
            "bladder ca": (load_hpo_by_name, "Bladder neoplasm"),
            "bowel cancer": (load_omim_by_id, 114500),
            "bowel polyps": (load_hpo_by_id, 200063),  # Colorectal polyposis
            "brain abnormalities": ABNORMAL_BRAIN,
            "brain abnormality": ABNORMAL_BRAIN,
            "brain malformation": ABNORMAL_BRAIN,
            "bulls ' eye maculopathy": BULLS_EYE_MACULOPATHY,  # TODO: Hacked due to us joining ' badly
            "caf au lait": (load_hpo_by_name, "Cafe-au-lait spot"),
            "carnitine transporter deficiency": (load_omim_by_id, 212140),
            "callosal dysgenesis": (load_hpo_by_name, 'Callosal agenesis'),
            "coagulation disorder": (load_hpo_by_name, "Abnormality of coagulation"),
            "congenital heart disease": (load_hpo_by_id, 1627),
            "congenital myasthenic": CMS,
            "congenital myasthenic syndrome": CMS,
            "congenital myaesthenic": CMS,
            "congenital myaesthenic syndrome": CMS,
            "cortical vision impairment": (load_hpo_by_id, 100704),  # Cerebral visual impairment
            "craniofacial dysmorphism": FACIAL_DYSMORPHISM,
            "crowded dentition": (load_hpo_by_id, 678),
            "development delay": DEVELOPMENTAL_DELAY,
            "dev issues": DEVELOPMENTAL_DELAY,
            "distal hypermobility": (load_hpo_by_name, "Limitation of joint mobility"),
            "duane syndrome": (load_hpo_by_id, 9921),
            "dystrophin": (load_gene_by_name, 'DMD'),
            "dysmorphic feature": DYSMORPHIC_FACE,
            "dysmorphic features": DYSMORPHIC_FACE,
            "easily bruised skin": (load_hpo_by_id, 978),  # Bruising susceptability
            "ehler danlos syndrome (type iii)": (load_omim_by_id, 130020),
            "ehlers-danos syndrome classic type": (load_omim_by_id, 130000),
            "elevated ammonia": (load_hpo_by_name, "Hyperammonemia"),
            "elevated ck": ELEVATED_CK,
            "elevated ketones": KETOSIS,
            "elevated lactate": (load_hpo_by_id, 2151),  # Increased serum lactate
            "elevated pth": (load_hpo_by_id, 3165),  # Elevated circulating parathyroid hormone
            "epileptic": (load_hpo_by_id, 1250),  # Seizures
            "facial dysmorphology": FACIAL_DYSMORPHISM,
            "fatty acid oxidation defect": FATTY_ACID_DISORDER,
            "fatty acid oxidation disorder": FATTY_ACID_DISORDER,
            "febrile sz": FEBRILE_SEIZURES,
            "fetal hydrops": HYDROPS_FETALIS,
            "global dd": DEVELOPMENTAL_DELAY,
            "global delay": DEVELOPMENTAL_DELAY,
            "global dev delay": DEVELOPMENTAL_DELAY,
            #"hailey-hailey syndrome" : (load_omim_by_name, "HAILEY-HAILEY DISEASE"),
            "hand flapping": (load_hpo_by_name, "Recurrent hand flapping"),
            "hearing aids": HEARING_IMPAIRMENT,
            "hearing impaired": HEARING_IMPAIRMENT,
            "hereditary neuralgic amyotrophy": (load_omim_by_id, 162100),
            "high ketones": KETOSIS,
            "high acth": (load_hpo_by_name, "Increased circulating ACTH level"),
            "hot flushes": (load_hpo_by_id, 32324),  # Episodic, so going for "Non-periodic recurrent fever"
            "hyperinsulinism": (load_hpo_by_id, 842),
            "hypoca": (load_hpo_by_name, "Hypocalcemia"),
            "hypoferritinaemia": (load_hpo_by_name, "Decreased serum ferritin"),  # hyper is there, hypo is not...
            "hypok": (load_hpo_by_name, "Hypokalemia"),
            "hypomg": (load_hpo_by_name, "Hypomagnesemia"),
            "hypop": (load_hpo_by_name, "Hypophosphatemia"),
            # I considered making a general conversion of "hypoplastic X" -> "Hypoplasia of X" but there are lots
            # of aliases that already do that, and a few entries for hypoplastic X but NOT hypoplasia of X so do case by case
            "hypoplastic right ventricle": (load_hpo_by_name, "Hypoplasia of right ventricle"),
            "inattention": (load_hpo_by_name, "Short attention span"),
            "increased renin": (load_hpo_by_id, 848),  # Increased circulating renin level
            "intellectual delay": (load_hpo_by_id, 1249),  # Intellectual disability (no delay anymore)
            "impaired consciousness": (load_hpo_by_name, "Reduced consciousness/confusion"),
            "iron deficiency": (load_hpo_by_id, 40130),  # Abnormal serum iron concentration
            "kneist dysplasia": (load_omim_by_name, "KNIEST DYSPLASIA"),
            "learning difficulties": COGNITIVE_IMPAIRMENT,
            "learning disability": COGNITIVE_IMPAIRMENT,
            "legius": (load_omim_by_name, "Legius Syndrome"),
            "leg pains": (load_hpo_by_name, "Limb pain"),
            "limb abnormalities": ABNORMALITY_OF_LIMBS,
            "low arylsulphatase": ARYLSULFATASE_A_DEFICIENCY,
            "low arylsulphatase A": ARYLSULFATASE_A_DEFICIENCY,
            "low bgl": (load_hpo_by_name, "Hypoglycemia"),
            "low bp": (load_hpo_by_id, 2615),  # Hypotension
            "low carnitine": (load_hpo_by_name, "Decreased plasma carnitine"),
            "lymphopaena": (load_hpo_by_name, "Lymphopenia"),
            "migranes": (load_hpo_by_name, "migraine"),
            "men type 1": (load_omim_by_id, 131100),
            "methylenetetrahyrofolate deficiency": (load_omim_by_id, 236250),  # HOMOCYSTINURIA DUE TO DEFICIENCY OF N(5,10)-METHYLENETETRAHYDROFOLATE REDUCTASE ACTIVITY
            "missing forearm": ABSENT_FOREARM,
            "missing forearms": ABSENT_FOREARM,
            "mitochondrial resp. chain disorder": MITO_DEFICIENCY,
            "mitochondrial respiratory chain disorder": MITO_DEFICIENCY,
            "moya moya": (load_hpo_by_id, 11834),
            "musculoskeletal abnormalities": (load_hpo_by_id, 33127),
            "na craving": (load_hpo_by_name, "Salt craving"),
            "neuroregression": (load_hpo_by_id, 2376),  # Developmental regression
            "noggin": (load_gene_by_name, 'NOG'),
            "no speech": (load_hpo_by_id, 1344),
            "ohtahara syndrome": (load_omim_by_id, 308350),
            "opisthoclonus": (load_hpo_by_name, "opisthotonus"),
            "opitz gbbb": (load_omim_by_id, 300000),
            "parkinson's disease": PARKINSONISM,
            "parkinsons": PARKINSONISM,
            "parkinson's": PARKINSONISM,
            "parkinson": PARKINSONISM,
            "pierre robin": PIERRE_ROBIN,
            "pierre-robin": PIERRE_ROBIN,
            'pig genes': PIG_GENES,
            "polysyndactyly": (load_hpo_by_id, 5873),
            "poor sleep": (load_hpo_by_id, 2360),
            "prolonged qt": (load_hpo_by_name, "Prolonged QT interval"),
            "prostate ca": (load_hpo_by_name, "Prostate cancer"),
            "pulmonary avms": PAVM,
            "pul avms": PAVM,
            "raised ck": ELEVATED_CK,
            "raised liver enzymes": (load_hpo_by_id, 2910),  # Elevated liver enzymes
            "raised ketones": KETOSIS,
            "raised methionine": (load_hpo_by_name, "Hypermethioninemia"),
            "raised urinary orotate": (load_hpo_by_name, "Oroticaciduria"),
            "raised tyrosine": (load_hpo_by_name, "Hypertyrosinemia"),
            "raised tsh": HIGH_TSH,
            "rem sleep": (load_hpo_by_id, 2494),  # Abnormal REM sleep
            "increased tsh": HIGH_TSH,
            "increased sweat": (load_hpo_by_name, "Hyperhidrosis"),
            "recurrent urtis": (load_hpo_by_name, "Recurrent upper respiratory tract infections"),
            "renal ca": (load_hpo_by_name, "Renal cell carcinoma"),
            "severe fetal hydrops": (load_hpo_by_name, "Severe hydrops fetalis"),
            "spastic cp": (load_hpo_by_name, "Cerebral palsy"),
            "thyroid ca": (load_hpo_by_name, "Thyroid carcinoma"),
            "type 1 diabetes": DIBETES_TYPE_1,
            "t1 diabetes": DIBETES_TYPE_1,
            "two hair whorls": (load_hpo_by_id, 10813),
            "uncoordinated": (load_hpo_by_id, 2406),
            "urea cycle": (load_genes_by_name, ["ARG1", "ASL", "ASS1", "CPS1", "NAGS", "OTC"]),
            "urogenital sinus": (load_hpo_by_id, 119),
            "waardenburg type ii": (load_omim_pks_containing_name, "waardenburg syndrome, type 2"),
            "widespread eyes": (load_hpo_by_id, 316),
        }

        # People put down eg Waardenburg but there are many different OMIM diseases - we'll put ALL of them
        # Switching to MONDO will help disease families, as it's hierarchial (unlike OMIM)
        ALPORT_SYNDROME = (load_omim_pks_containing_name, "Alport Syndrome")
        CILIARY_DYSKINESIA = (load_omim_pks_containing_name, "Ciliary dyskinesia")
        BARTTER_SYNDROME = (load_omim_pks_containing_name, "Bartter Syndrome")
        BRUGADA_SYNDROME = (load_omim_pks_containing_name, "Brugada Syndrome")
        CHARCOT_MARIE_TOOTH = (load_omim_pks_containing_name, "Charcot-marie-tooth")
        CHONDRODYSPLASIA_PUNCTATA = (load_omim_pks_containing_name, "CHONDRODYSPLASIA PUNCTATA")
        EHLER_DANOS = (load_omim_pks_containing_name, "EHLERS-DANLOS SYNDROME")
        HLH = (load_omim_pks_containing_name, "HEMOPHAGOCYTIC LYMPHOHISTIOCYTOSIS")
        LEBER_AMAUROSIS = (load_omim_pks_containing_name, "leber congenital amaurosis")
        ROBINOW = (load_omim_pks_containing_name, "Robinow Syndrome")
        SANFILIPPO = (load_omim_pks_containing_alias_name, "Sanfilippo Syndrome")
        STICKLER = (load_omim_pks_containing_name, "Stickler Syndrome")
        TUBEROUS_SCLEROSIS = (load_omim_pks_containing_name, "tuberous sclerosis")

        DISEASE_FAMILIES = {
            "aicardi goutires syndrome": (load_omim_pks_containing_name, "AICARDI-GOUTIERES"),
            "alport syndrome": ALPORT_SYNDROME,
            "alport": ALPORT_SYNDROME,
            "autoinflammatory syndrome": (load_omim_pks_containing_name, "autoinflammatory syndrome"),
            "bartter syndrome": BARTTER_SYNDROME,
            "bartter": BARTTER_SYNDROME,
            "brugada": BRUGADA_SYNDROME,
            "brugada syndrome": BRUGADA_SYNDROME,
            "charcot-marie-tooth": CHARCOT_MARIE_TOOTH,
            "charcot marie tooth": CHARCOT_MARIE_TOOTH,
            "chondrodysplasia": CHONDRODYSPLASIA_PUNCTATA,
            "chondroplasia punctata": CHONDRODYSPLASIA_PUNCTATA,  # Typo
            "cilial dyskinesis": CILIARY_DYSKINESIA,
            "cilial dyskinesia": CILIARY_DYSKINESIA,
            "ehrlers danlos": EHLER_DANOS,
            "ehlers-danos": EHLER_DANOS,
            "ehler danlos": EHLER_DANOS,
            "gaucher disease":  (load_omim_pks_containing_name, "GAUCHER DISEASE"),
            "glycogen storage disease": GLYCOGEN_STORAGE_DISEASE,
            "glut1 deficiency": (load_omim_pks_containing_name, "GLUT1 DEFICIENCY SYNDROME"),
            "hemophagocytic lymphohistiocytosis": HLH,
            "hlh": HLH,
            "hht": (load_omim_pks_containing_name, "Hereditary hemorrhagic telangiectasia"),
            "leber amaurosis": LEBER_AMAUROSIS,
            "lebers amaurosis": LEBER_AMAUROSIS,
            "leber's amaurosis": LEBER_AMAUROSIS,  # TODO: This doesn't match as seems ' ' inserted??
            "osteogenesis imperfecta": (load_omim_pks_containing_name, "osteogenesis imperfecta"),
            "robinow syndrome": ROBINOW,
            "robinow": ROBINOW,
            "sanfilippo": SANFILIPPO,
            "stickler syndrome": STICKLER,
            "stickler": STICKLER,
            "tuberous sclerosis": TUBEROUS_SCLEROSIS,
            "usher syndrome": (load_omim_pks_containing_name, "usher syndrome"),
            "waardenburg": (load_omim_pks_containing_name, "waardenburg"),
        }
        return HARDCODED_LOOKUPS, CASE_INSENSITIVE_LOOKUPS, DISEASE_FAMILIES
