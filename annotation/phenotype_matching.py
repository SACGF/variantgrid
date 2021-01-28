from collections import defaultdict
from typing import List, Optional, Dict, Union, Tuple, Iterable, Any

from django.db.models import Count
from django.db.models.functions import Lower
import Levenshtein
import logging
import nltk
import re
import time


from annotation.models.models_phenotype_match import PhenotypeMatchTypes, \
    TextPhenotypeMatch, PhenotypeDescription, TextPhenotype, TextPhenotypeSentence
from genes.models import GeneSymbol, Gene
from library.log_utils import log_traceback
from library.utils import get_and_log_time_since, invert_dict_of_lists
from ontology.models import OntologyTerm, OntologyService
from patients.models import Patient

HPO_PATTERN = re.compile(r"HP:(\d{7})$")
HPO_TYPO_PATTERN = re.compile(r"HPO:(\d{7})$")
OMIM_PATTERN = re.compile(r"OMIM:(\d+)$")

MIN_MATCH_LENGTH = 3
MIN_LENGTH_SINGLE_WORD_FUZZY_MATCH = 5
MAX_COMBO_LENGTH = 14  # Checked HPO words in DB

CodePK = Any
Lookups = Dict[CodePK, str]
OntologyObj = Union[OntologyTerm, Gene]
OntologyDict = Dict[str, OntologyObj]
OntologyResults = Tuple[str, List[Union[CodePK, OntologyObj]]]


def get_word_combos_and_spans(words_and_spans: List, max_combo_length: Optional[int]) -> List:
    combos_and_spans = []
    for i in range(len(words_and_spans)):
        combo_length = len(words_and_spans) - i
        if max_combo_length:
            combo_length = min(max_combo_length, combo_length)
        for j in range(i, i + combo_length):
            combo_and_span = words_and_spans[i:j + 1]
            combos_and_spans.append(combo_and_span)

    return combos_and_spans


def get_word_combos_and_spans_sorted_by_length(words_and_spans, max_combo_length: Optional[int] = None) -> Iterable:
    word_combos_and_spans = get_word_combos_and_spans(words_and_spans, max_combo_length)
    return reversed(sorted(word_combos_and_spans, key=lambda item: sum([len(i[0]) for i in item])))


def get_id_from_multi_word_fuzzy_match(lookup: Lookups, words: List[str], text: str, distance: int) -> Optional[CodePK]:
    potentials: Lookups = dict()
    for w in words:
        potentials.update(lookup[w])

    if not potentials:
        return None
    return get_id_from_fuzzy_match(potentials, text, distance)


def get_id_from_single_word_fuzzy_match(single_words_by_length: Dict[int, Lookups], text: str, distance: int) -> Optional[CodePK]:
    text_length = len(text)
    potentials: Lookups = dict()
    for l in [text_length - distance, text_length, text_length + distance]:
        words = single_words_by_length.get(l)
        if words:
            potentials.update(words)
    return get_id_from_fuzzy_match(potentials, text, distance)


def get_id_from_fuzzy_match(lookup: Lookups, text: str, max_distance: int) -> Optional[CodePK]:
    for description, pk in lookup.items():
        distance = Levenshtein.distance(description, text)  # @UndefinedVariable
        #print("'%s' <-> '%s' distance: %d" % (description, text, distance))
        if distance <= max_distance:
            return pk
    return None


def get_multi_word_hpo_fuzzy(hpo_word_lookup, words, text: str, distance: int = 1) -> Optional[CodePK]:
    return get_id_from_multi_word_fuzzy_match(hpo_word_lookup, words, text, distance)


def get_multi_word_omim_fuzzy(omim_word_lookup, words, text: str, distance=1) -> Optional[CodePK]:
    return get_id_from_multi_word_fuzzy_match(omim_word_lookup, words, text, distance)


def get_single_word_hpo_fuzzy(hpo_single_words_by_length: Dict[int, Lookups], text: str, distance=1) -> Optional[CodePK]:
    return get_id_from_single_word_fuzzy_match(hpo_single_words_by_length, text, distance)


def get_single_word_omim_fuzzy(omim_single_words_by_length: Dict[int, Lookups], text: str, distance=1) -> Optional[CodePK]:
    return get_id_from_single_word_fuzzy_match(omim_single_words_by_length, text, distance)


def create_word_lookups(records: Lookups) -> Dict[str, OntologyDict]:
    word_lookup = defaultdict(dict)

    for text, obj in records.items():
        for word in text.split():
            word_lookup[word][text] = obj

    return word_lookup


def get_special_case_match(text, hpo_records: OntologyDict, omim_records: OntologyDict, gene_records: OntologyDict) -> Tuple[List[OntologyObj], List[OntologyObj], List[OntologyObj]]:
    def load_omim_alias_by_id(accession) -> OntologyResults:
        pk = OntologyService.index_to_id(OntologyService.OMIM, accession)
        return PhenotypeMatchTypes.OMIM, [pk]

    def load_omim_by_name(description: str) -> OntologyResults:
        omim_alias = omim_records[description.lower()]
        return PhenotypeMatchTypes.OMIM, [omim_alias]

    def load_hpo_by_id(hpo_id) -> OntologyResults:
        pk = OntologyService.index_to_id(OntologyService.HPO, hpo_id)
        hpo = OntologyTerm.objects.get(pk=pk)
        return PhenotypeMatchTypes.HPO, [hpo.pk]

    def load_hpo_by_name(hpo_name) -> OntologyResults:
        hpo = hpo_records[hpo_name.lower()]
        return PhenotypeMatchTypes.HPO, [hpo]

    def load_hpo_list_by_names(hpo_name_list) -> OntologyResults:
        hpo_list: List[OntologyObj] = list()
        for hpo_name in hpo_name_list:
            _, hpo = load_hpo_by_name(hpo_name)
            hpo_list.extend(hpo)
        return PhenotypeMatchTypes.HPO, hpo_list

    def load_gene_by_name(gene_symbol: str) -> OntologyResults:
        gene = gene_records[gene_symbol.lower()]
        return PhenotypeMatchTypes.GENE, [gene]

    def load_genes_by_name(gene_symbols_list: List[str]) -> OntologyResults:
        genes_list = []
        for gene_symbol in gene_symbols_list:
            _, genes = load_gene_by_name(gene_symbol)
            genes_list.extend(genes)

        return PhenotypeMatchTypes.GENE, genes_list

    omim_qs = OntologyTerm.objects.filter(ontology_service=OntologyService.OMIM)

    def load_omim_pks_containing_name(name: str) -> OntologyResults:
        return PhenotypeMatchTypes.OMIM, omim_qs.filter(name__icontains=name).values_list("pk", flat=True)

    def load_omim_pks_containing_alias_name(name) -> OntologyResults:
        return PhenotypeMatchTypes.OMIM, omim_qs.objects.filter(aliases__icontains=name).values_list("pk", flat=True)

    ABSENT_FOREARM = (load_hpo_by_name, 'absent forearm')
    ABNORMAL_BRAIN = (load_hpo_by_name, "Abnormality of brain morphology")
    ABNORMALITY_OF_LIMBS = (load_hpo_by_name, 'Abnormality of limbs')
    ARYLSULFATASE_A_DEFICIENCY = (load_omim_by_name, "ARYLSULFATASE A DEFICIENCY")
    AUTISTIC = (load_hpo_by_id, 729)
    BULLS_EYE_MACULOPATHY = (load_hpo_by_name, "bull's eye maculopathy")
    FATTY_ACID_DISORDER = (load_hpo_by_name, "Abnormality of fatty-acid metabolism")
    HUS = (load_hpo_by_name, "Hemolytic-uremic syndrome")
    MITO_DEFICIENCY = (load_omim_by_name, "MITOCHONDRIAL COMPLEX I DEFICIENCY")
    DEVELOPMENTAL_DELAY = (load_hpo_by_name, "Developmental delay")
    GLOBAL_DEVELOPMENTAL_DELAY = (load_hpo_by_name, "Global developmental delay")
    ELEVATED_CK = (load_hpo_by_name, "Elevated creatine kinase")
    KETOSIS = (load_hpo_by_name, "Ketosis")
    PIERRE_ROBIN = (load_hpo_by_name, "Pierre-Robin sequence")
    PAVM = (load_hpo_by_name, "Pulmonary arteriovenous malformation")
    CMS = (load_hpo_by_name, "Fatigable weakness")
    HYDROPS_FETALIS = (load_hpo_by_name, "Nonimmune hydrops fetalis")
    AFEBRILE = (load_hpo_by_name, "Focal seizures, afebril")
    GEFS = (load_hpo_by_name, "Febrile seizures")  # GEFS+ is a multi-type OMIM disease, this links to all those though
    PIG_GENES = (load_genes_by_name, ['PIG' + i for i in 'ABCFGHKLMNOPQSTUVWXYZ'])
    GLYCOGEN_STORAGE_DISEASE = (load_omim_pks_containing_name, "glycogen storage disease")
    HIGH_TSH = (load_hpo_by_name, 'Thyroid-stimulating hormone excess')
    PARKINSONISM = (load_hpo_by_name, 'Parkinsonism')
    DIBETES_TYPE_1 = (load_hpo_by_name, 'Type I diabetes mellitus')
    DYSMORPHIC_FACE = (load_hpo_by_name, "Abnormal facial shape")
    COGNITIVE_IMPAIRMENT = (load_hpo_by_name, "Cognitive impairment")  # There is also 'Specific learning disability' but this is different I think
    FACIAL_DYSMORPHISM = (load_hpo_by_id, 1999)
    HEARING_IMPAIRMENT = (load_hpo_by_name, "Hearing impairment")

    HARDCODED_LOOKUPS = {'aHUS': HUS,
                         "ALL": (load_hpo_by_name, "Acute lymphoblastic leukemia"),
                         # AML fix until we get new HPO data - see https://github.com/obophenotype/human-phenotype-ontology/issues/4236
                         "AML": (load_hpo_by_name, "Acute myeloid leukemia"),
                         "ADPCKD": (load_omim_by_name, "POLYCYSTIC KIDNEY DISEASE 1"),
                         "AVSD": (load_hpo_by_name, "Atrioventricular septal defect"),
                         "BCC": (load_hpo_by_name, "Basal cell carcinoma"),
                         "BrCa": (load_omim_alias_by_id, 114480),  # BREAST CANCER
                         "CHD": (load_hpo_by_name, "Abnormal heart morphology"),
                         "CMS": CMS,
                         "DD":  DEVELOPMENTAL_DELAY,
                         "FAOD": (load_hpo_by_name, "Abnormality of fatty-acid metabolism"),
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
                         "IBD": (load_omim_alias_by_id, 266600),  # IBD1
                         "ID": (load_hpo_by_name, 'intellectual disability'),
                         "LGA": (load_hpo_by_name, "Large for gestational age"),
                         "LQTS": (load_hpo_by_name, "Long QT syndrome"),
                         "MM": (load_hpo_by_name, 'Multiple myeloma'),
                         "NCS": (load_hpo_by_name, "Neurocardiogenic syncope"),
                         "PCKD": (load_hpo_by_name, "Polycystic kidney dysplasia"),
                         "PV": (load_omim_alias_by_id, 263300),  # POLYCYTHEMIA VERA; PV
                         "SCID": (load_hpo_by_name, "Severe combined immunodeficiency"),
                         'SMA': (load_hpo_by_name, "spinal muscular atrophy"),
                         "SNA12": (load_gene_by_name, "SNAI2"),  # Common misspelling
                         "SUDEP": (load_hpo_by_name, ["Sudden death", "Epilepsy"]),
                         "VSD": (load_hpo_by_name, "Ventricular septal defect")}

    CASE_INSENSITIVE_LOOKUPS = {"aarskog": (load_omim_by_name, "AARSKOG-SCOTT SYNDROME"),
                                "abdo pain": (load_hpo_by_name, "Abdominal pain"),
                                "abnormal mri brain": ABNORMAL_BRAIN,
                                "aching limbs": (load_hpo_by_name, "Limb pain"),
                                "adenosine phosphoribosyl transferase deficiencies": (load_omim_alias_by_id, 614723),
                                "agenesis cc": (load_hpo_by_name, "Agenesis of corpus callosum"),
                                "afebrile seizures": AFEBRILE,
                                "afebrile": AFEBRILE,
                                "autistic features": AUTISTIC,
                                "autistic": AUTISTIC,
                                "behaviour problems": (load_hpo_by_name, "Behavioral abnormality"),
                                "bladder ca": (load_hpo_by_name, "Bladder neoplasm"),
                                "bowel cancer": (load_omim_alias_by_id, 114500),
                                "bowel polyps": (load_hpo_by_name, "Colorectal polyps"),
                                "brain abnormalities": ABNORMAL_BRAIN,
                                "brain abnormality": ABNORMAL_BRAIN,
                                "brain malformation": ABNORMAL_BRAIN,
                                "bulls ' eye maculopathy": BULLS_EYE_MACULOPATHY,  # TODO: Hacked due to us joining ' badly
                                "caf au lait": (load_hpo_by_name, "Cafe-au-lait spot"),
                                "carnitine transporter deficiency": (load_omim_alias_by_id, 212140),
                                "callosal dysgenesis": (load_hpo_by_name, 'Callosal agenesis'),
                                "coagulation disorder": (load_hpo_by_name, "Abnormality of coagulation"),
                                "congenital heart disease": (load_hpo_by_id, 1627),
                                "congenital myasthenic": CMS,
                                "congenital myasthenic syndrome": CMS,
                                "congenital myaesthenic": CMS,
                                "congenital myaesthenic syndrome": CMS,
                                "cortical vision impairment": (load_hpo_by_name, "Cortical visual impairment"),
                                "craniofacial dysmorphism": FACIAL_DYSMORPHISM,
                                "crowded dentition": (load_hpo_by_id, 678),
                                "development delay": DEVELOPMENTAL_DELAY,
                                "dev issues": DEVELOPMENTAL_DELAY,
                                "distal hypermobility": (load_hpo_by_name, "Limitation of joint mobility"),
                                "duane syndrome": (load_hpo_by_id, 9921),
                                "dystrophin": (load_gene_by_name, 'DMD'),
                                "dysmorphic feature": DYSMORPHIC_FACE,
                                "dysmorphic features": DYSMORPHIC_FACE,
                                "easily bruised skin": (load_hpo_by_name, "Bruise easily"),
                                "ehler danlos syndrome (type iii)": (load_omim_alias_by_id, 130020),
                                "ehlers-danos syndrome classic type": (load_omim_alias_by_id, 130000),
                                "elevated ammonia": (load_hpo_by_name, "Hyperammonemia"),
                                "elevated ck": ELEVATED_CK,
                                "elevated ketones": KETOSIS,
                                "elevated lactate": (load_hpo_by_name, "Increased blood lactate"),
                                "elevated pth": (load_hpo_by_name, "Elevated circulating parathyroid hormone (PTH) level"),
                                "epileptic": (load_hpo_by_name, "epilepsy"),
                                "facial dysmorphology": FACIAL_DYSMORPHISM,
                                "fatty acid oxidation defect": FATTY_ACID_DISORDER,
                                "fatty acid oxidation disorder": FATTY_ACID_DISORDER,
                                "fetal hydrops": HYDROPS_FETALIS,
                                "global dd": GLOBAL_DEVELOPMENTAL_DELAY,
                                "global delay": GLOBAL_DEVELOPMENTAL_DELAY,
                                "global dev delay": GLOBAL_DEVELOPMENTAL_DELAY,
                                #"hailey-hailey syndrome" : (load_omim_by_name, "HAILEY-HAILEY DISEASE"),
                                "hand flapping": (load_hpo_by_name, "Recurrent hand flapping"),
                                "hearing aids": HEARING_IMPAIRMENT,
                                "hearing impaired": HEARING_IMPAIRMENT,
                                "hereditary neuralgic amyotrophy": (load_omim_by_name, "AMYOTROPHY, HEREDITARY NEURALGIC"),
                                "high ketones": KETOSIS,
                                "high acth": (load_hpo_by_name, "Increased circulating ACTH level"),
                                "hot flushes": (load_hpo_by_name, "Episodic fever"),  # Not the same but best I can match
                                "hyperinsulinism": (load_hpo_by_name, "Elevated insulin level"),
                                "hypoca": (load_hpo_by_name, "Hypocalcemia"),
                                "hypoferritinaemia": (load_hpo_by_name, "Decreased serum ferritin"),  # hyper is there, hypo is not...
                                "hypok": (load_hpo_by_name, "Hypokalemia"),
                                "hypomg": (load_hpo_by_name, "Hypomagnesemia"),
                                "hypop": (load_hpo_by_name, "Hypophosphatemia"),
                                # I considered making a general conversion of "hypoplastic X" -> "Hypoplasia of X" but there are lots
                                # of aliases that already do that, and a few entries for hypoplastic X but NOT hypoplasia of X so do case by case
                                "hypoplastic right ventricle": (load_hpo_by_name, "Hypoplasia of right ventricle"),
                                "inattention": (load_hpo_by_name, "Short attention span"),
                                "increased renin": (load_hpo_by_name, "Increased serum renin"),
                                "intellectual delay": (load_hpo_by_name, "Delayed intellectual development"),
                                "impaired consciousness": (load_hpo_by_name, "Reduced consciousness/confusion"),
                                "iron deficiency": (load_hpo_by_name, "Abnormal serum iron"),
                                "kneist dysplasia": (load_omim_by_name, "KNIEST DYSPLASIA"),
                                "learning difficulties": COGNITIVE_IMPAIRMENT,
                                "learning disability": COGNITIVE_IMPAIRMENT,
                                "legius": (load_omim_by_name, "Legius Syndrome"),
                                "leg pains": (load_hpo_by_name, "Limb pain"),
                                "limb abnormalities": ABNORMALITY_OF_LIMBS,
                                "low arylsulphatase": ARYLSULFATASE_A_DEFICIENCY,
                                "low arylsulphatase A": ARYLSULFATASE_A_DEFICIENCY,
                                "low bgl": (load_hpo_by_name, "Hypoglycemia"),
                                "low bp": (load_hpo_by_name, "Low blood pressure"),
                                "low carnitine": (load_hpo_by_name, "Decreased plasma carnitine"),
                                "lymphopaena": (load_hpo_by_name, "Lymphopenia"),
                                "migranes": (load_hpo_by_name, "migraine"),
                                "men type 1": (load_omim_alias_by_id, 131100),
                                "methylenetetrahyrofolate deficiency": (load_omim_alias_by_id, 236250),  # HOMOCYSTINURIA DUE TO DEFICIENCY OF N(5,10)-METHYLENETETRAHYDROFOLATE REDUCTASE ACTIVITY
                                "missing forearm": ABSENT_FOREARM,
                                "missing forearms": ABSENT_FOREARM,
                                "mitochondrial resp. chain disorder": MITO_DEFICIENCY,
                                "mitochondrial respiratory chain disorder": MITO_DEFICIENCY,
                                "moya moya": (load_hpo_by_id, 11834),
                                "musculoskeletal abnormalities": (load_hpo_list_by_names, ["Muscular abnormality", "Skeletal abnormalities"]),
                                "na craving": (load_hpo_by_name, "Salt craving"),
                                "neuroregression": (load_hpo_by_name, "Neurodevelopmental regression"),
                                "noggin": (load_gene_by_name, 'NOG'),
                                "no speech": (load_hpo_by_id, 1344),
                                "ohtahara syndrome": (load_omim_by_name, "OHTAHARA SYNDROME, X-LINKED"),
                                "opisthoclonus": (load_hpo_by_name, "opisthotonus"),
                                "opitz gbbb": (load_omim_by_name, "OPITZ GBBB SYNDROME, X-LINKED"),
                                "parkinson's disease": PARKINSONISM,
                                "parkinsons": PARKINSONISM,
                                "parkinson's": PARKINSONISM,
                                "parkinson": PARKINSONISM,
                                "pierre robin": PIERRE_ROBIN,
                                "pierre-robin": PIERRE_ROBIN,
                                'pig genes': PIG_GENES,
                                "periodic fever": (load_omim_by_name, "PERIODIC FEVER, FAMILIAL, AUTOSOMAL DOMINANT"),
                                "polysyndactyly": (load_hpo_by_name, "Polysyndactyly of big toe"),
                                "poor sleep": (load_hpo_by_id, 2360),
                                "prolonged qt": (load_hpo_by_name, "Prolonged QT interval"),
                                "prostate ca": (load_hpo_by_name, "Prostate cancer"),
                                "pulmonary avms": PAVM,
                                "pul avms": PAVM,
                                "raised ck": ELEVATED_CK,
                                "raised liver enzymes": (load_hpo_by_name, "Elevated liver enzymes"),
                                "raised ketones": KETOSIS,
                                "raised methionine": (load_hpo_by_name, "Hypermethioninemia"),
                                "raised urinary orotate": (load_hpo_by_name, "High urine orotic acid levels"),
                                "raised tyrosine": (load_hpo_by_name, "Hypertyrosinemia"),
                                "raised tsh": HIGH_TSH,
                                "increased tsh": HIGH_TSH,
                                "increased sweat": (load_hpo_by_name, "Hyperhidrosis"),
                                "recurrent urtis": (load_hpo_by_name, "Recurrent upper respiratory tract infections"),
                                "rem sleep": (load_hpo_by_name, "Abnormal REM sleep"),
                                "renal ca": (load_hpo_by_name, "Renal cell carcinoma"),
                                "severe fetal hydrops": (load_hpo_by_name, "Severe hydrops fetalis"),
                                "spastic cp": (load_hpo_by_name, "Cerebral palsy"),
                                "thyroid ca": (load_hpo_by_name, "Thyroid carcinoma"),
                                "type 1 diabetes": DIBETES_TYPE_1,
                                "t1 diabetes": DIBETES_TYPE_1,
                                "two hair whorls": (load_hpo_by_id, 10813),
                                "uncoordinated": (load_hpo_by_id, 2406),
                                "urea cycle": (load_genes_by_name, ["ARG1", "ASL", "ASS1", "CPS1", "NAGS", "OTC"]),
                                "urogenital sinus": (load_hpo_by_name, 'Urogenital anomalies'),
                                "waardenburg type ii": (load_omim_pks_containing_name, "waardenburg syndrome, type 2"),
                                "widespread eyes": (load_hpo_by_name, "Widely spaced eyes")}

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

    hpo_list = []
    omim_alias_list = []
    genes = []

    hl = HARDCODED_LOOKUPS.get(text)
    if not hl:
        # Lookup accessions in form of: HP:0000362, OMIM:607196

        PATTERNS = [
            (HPO_PATTERN, load_hpo_by_id),
            (HPO_TYPO_PATTERN, load_hpo_by_id),
            (OMIM_PATTERN, load_omim_alias_by_id),
        ]

        for pattern, load_func in PATTERNS:
            if m := pattern.match(text):
                accession = int(m.group(1))
                hl = (load_func, accession)
                break

    if not hl:  # Try lowercase lookups
        lowercase_text = text.lower()

        hl = CASE_INSENSITIVE_LOOKUPS.get(lowercase_text)
        if not hl:
            hl = DISEASE_FAMILIES.get(lowercase_text)

    if hl:
        (func, arg) = hl
        try:
            match_type, records = func(arg)
            if match_type == PhenotypeMatchTypes.HPO:
                hpo_list.extend(records)
            elif match_type == PhenotypeMatchTypes.OMIM:
                omim_alias_list.extend(records)
            elif match_type == PhenotypeMatchTypes.GENE:
                genes.extend(records)

            #logging.info("Got exact: %s => %s", text, records)
        except Exception as e:
            msg = f"Error: {e}, func: {func}, arg={arg}"
            logging.error(msg)
            log_traceback()
#            raise ValueError(msg)

    return hpo_list, omim_alias_list, genes


class SkipAllPhenotypeMatchException(Exception):
    pass


def words_together(text, first_words, second_words):
    first_words = {x.lower() for x in first_words}
    second_words = {x.lower() for x in second_words}

    for f in first_words:
        if f in text:
            for s in second_words:
                if s in text:
                    return True
    return False


def skip_word(lower_text):
    """ Return true to skip a word, throws SkipAllPhenotypeMatchException to skip all. Only need to skip >MIN_LENGTH words as will do that later (after exact) """

    # For multi-words where you want to skip components
    SKIP_ALL = {"library prep", "to cgf", "set up", "ad pattern", "recurrent eps", "rest of"}
    if lower_text in SKIP_ALL:
        raise SkipAllPhenotypeMatchException()

    if words_together(lower_text, {"TAT"}, {"non-urgent", "months", "month", "days", "week", "weeks"}):
        raise SkipAllPhenotypeMatchException()

    if words_together(lower_text, {"trio"}, {"exome", "MedEx", "WES", "TS1", "father", "mother"}):
        raise SkipAllPhenotypeMatchException()

    # Words which have no use matching on their own
    COMMON_WORDS = {'acute', 'adult', 'all', 'and', 'associated',
                    'bad', 'bilateral', 'birth', 'blood', 'borderline',
                    'can', 'carries', 'central', 'change', 'charge', 'child', 'chronic', 'close', 'comma', 'commas', 'coned', 'cord', 'cousin', 'cousins',
                    'diffused', 'deficiency', 'disorder', 'distal',
                    'exclude',
                    'face', 'familial', 'floating', 'focal', 'forms', 'frequent', 'frequency', 'from',
                    'generalized', "generalised",
                    'hard', 'hearing',
                    'image', 'inheritance', 'insulin',
                    'joints',
                    'kit',
                    'large', 'lateral', 'left', 'likes', 'liver',
                    'match', 'macro', 'mild', 'milena', 'moderate', 'motor',
                    'name',
                    'onset',
                    'parts', 'pending', 'periodic', 'person', 'pit', 'plan', 'position', 'profound', 'prolonged', 'proximal', 'progressive',
                    'range', 'raise', 'recurrent', 'right',
                    'severe', 'she', 'short', 'skeletal', 'sleep', 'syndrome',
                    'the', 'transient',
                    'wants', 'was', 'with'}

    return lower_text in COMMON_WORDS


def calculate_match_distance(words: List[str]) -> int:
    """ by default we match on 1 - however we may want to be a bit lax sometimes """
    num_ae_words = 0
    for w in words:
        num_ae_words += "ae" in w
    distance = max(1, num_ae_words)
    return distance


def get_terms_from_words(
        text_phenotype,
        words_and_spans_subset,
        hpo_records,
        hpo_word_lookup,
        hpo_single_words_by_length,
        omim_records,
        omim_word_lookup,
        omim_single_words_by_length,
        gene_records):
    words = [ws[0] for ws in words_and_spans_subset]
    text = ' '.join(words)
    lower_text = text.lower()

    hpo_list = []
    omim_list = []
    gene_symbols = []
    special_case_match = get_special_case_match(text, hpo_records, omim_records, gene_records)
    if any(special_case_match):
        hpo_list, omim_list, gene_symbols = special_case_match
    else:
        if len(lower_text) < MIN_MATCH_LENGTH:
            return []

        if skip_word(lower_text):
            return []

        hpo = hpo_records.get(lower_text)
        omim = omim_records.get(lower_text)
        if not any([hpo, omim]):
            if len(words) == 1:
                w = words[0]
                if len(w) >= MIN_LENGTH_SINGLE_WORD_FUZZY_MATCH:
                    hpo = get_single_word_hpo_fuzzy(hpo_single_words_by_length, lower_text)
                    omim = get_single_word_omim_fuzzy(omim_single_words_by_length, lower_text)
            else:
                lower_words = [w.lower() for w in words]
                distance = calculate_match_distance(lower_words)
                hpo = get_multi_word_hpo_fuzzy(hpo_word_lookup, lower_words, lower_text, distance=distance)
                omim = get_multi_word_omim_fuzzy(omim_word_lookup, lower_words, lower_text, distance=distance)

        if len(words) == 1:
            # Don't do fuzzy for genes as likely to get false positives
            gene_symbol = gene_records.get(lower_text)
            if gene_symbol:
                gene_symbols.append(gene_symbol)

        if hpo:
            hpo_list.append(hpo)

        if omim:
            omim_list.append(omim)

    results = []
    offset_start = words_and_spans_subset[0][1][0]
    offset_end = words_and_spans_subset[-1][1][1]

    for ontology_term_id in hpo_list:
        tpm = TextPhenotypeMatch.objects.create(text_phenotype=text_phenotype,
                                                match_type=PhenotypeMatchTypes.HPO,
                                                ontology_term_id=ontology_term_id,
                                                offset_start=offset_start,
                                                offset_end=offset_end)
        results.append(tpm)

    for ontology_term_id in omim_list:
        tpm = TextPhenotypeMatch.objects.create(text_phenotype=text_phenotype,
                                                match_type=PhenotypeMatchTypes.OMIM,
                                                ontology_term_id=ontology_term_id,
                                                offset_start=offset_start,
                                                offset_end=offset_end)
        results.append(tpm)

    for gene_symbol_id in gene_symbols:
        tpm = TextPhenotypeMatch.objects.create(text_phenotype=text_phenotype,
                                                match_type=PhenotypeMatchTypes.GENE,
                                                gene_symbol_id=gene_symbol_id,
                                                offset_start=offset_start,
                                                offset_end=offset_end)
        results.append(tpm)

    return results


def sub_array_index(array, sub_array):
    for i in range(0, len(array) - len(sub_array) + 1):
        if array[i:i + len(sub_array)] == sub_array:
            return i
    return None


def parse_words(text_phenotype,
                input_words_and_spans,
                hpo_records,
                hpo_word_lookup,
                hpo_single_words_by_length,
                omim_records,
                omim_word_lookup,
                omim_single_words_by_length,
                gene_records) -> List:
    results = []

    word_combos_and_spans = get_word_combos_and_spans_sorted_by_length(input_words_and_spans, max_combo_length=MAX_COMBO_LENGTH)
    for words_and_spans_subset in word_combos_and_spans:
        words_results = get_terms_from_words(text_phenotype, words_and_spans_subset, hpo_records, hpo_word_lookup, hpo_single_words_by_length, omim_records, omim_word_lookup, omim_single_words_by_length, gene_records)
        if words_results:
            results.extend(words_results)
            if words_and_spans_subset != input_words_and_spans:  # More to match
                i = sub_array_index(input_words_and_spans, words_and_spans_subset)
                before_words = input_words_and_spans[:i]
                after_words = input_words_and_spans[i + len(words_and_spans_subset):]

                if before_words:
                    before_results = parse_words(text_phenotype, before_words, hpo_records, hpo_word_lookup, hpo_single_words_by_length, omim_records, omim_word_lookup, omim_single_words_by_length, gene_records)
                    results.extend(before_results)

                if after_words:
                    after_results = parse_words(text_phenotype, after_words, hpo_records, hpo_word_lookup, hpo_single_words_by_length, omim_records, omim_word_lookup, omim_single_words_by_length, gene_records)
                    results.extend(after_results)

            break

    return results


def split_adj_noun_and_noun(words, tags, spans):
    """ This is used to transform eg "long fingers and toes" into "long fingers and long toes"
        I looked and this pattern doesn't appear in HPO terms so is ok to split up """

    new_words_and_spans = []
    if tags == ['JJ', 'NNS', 'CC', 'NNS'] and words[2].lower() == 'and':
        new_order = [0, 1, 2, 0, 3]
        for i in new_order:
            new_words_and_spans.append((words[i], spans[i]))

    return new_words_and_spans


def transform_words(words, tags, spans):
    # TODO:
    # X disorder => Abnormality of X
    # X Ca => X Cancer?

    TRANSFORMATIONS = [split_adj_noun_and_noun]

    for t in TRANSFORMATIONS:
        new_words_and_spans = t(words, tags, spans)
        if new_words_and_spans:
            return new_words_and_spans

    return list(zip(words, spans))


def tokens_and_spans(txt, tokens):
    offset = 0

    tokenized_text_and_spans = []
    for token in tokens:
        offset = txt.find(token, offset)
        tokenized_text_and_spans.append((token, offset, offset + len(token)))
        offset += len(token)

    return tokenized_text_and_spans


def handle_unsplit(unsplit_digits_and_spans, split_words_and_spans):
    if unsplit_digits_and_spans:
        unsplit_word = '/'.join([d[0] for d in unsplit_digits_and_spans])
        my_offset_start = unsplit_digits_and_spans[0][1]
        my_offset_end = unsplit_digits_and_spans[-1][2]
        split_words_and_spans.append((unsplit_word, my_offset_start, my_offset_end))
        unsplit_digits_and_spans[:] = []


def initial_replace_characters(text):
    REPLACE_CHARS = {';': ' '}

    cleaned_text = text
    for from_char, to_char in REPLACE_CHARS.items():
        cleaned_text = cleaned_text.replace(from_char, to_char)

    return cleaned_text


def word_tokenise_and_spans(txt):
    """ Split by nltk.sent_tokenize and / (only if alphanumeric not for dates) """

    cleaned_text = initial_replace_characters(txt)
    tokens = nltk.word_tokenize(cleaned_text)
    split_words_and_spans = []

    # {TO: FROM} - These have to be the same size! Then switch them out later
    PRE_SPLIT_REPLACE = {"dv/p": "dv_p"}

    # replace word : List of alternate words to replace
    REPLACE_ALTERNATES = {
        "bartter": ["bartters", "bartter's"],
        "bull's eye": ["bullseye"],
        "developmental": ["dev", "devt", "dv_p", "devel"],
        "dysmorphic": ["dysm"],
        "gaucher": ["gauchers", "gaucher's"],
        "hemophagocytic": ["haemophagcytic", "haemophogocytic", "haemophagocytic"],
        "leukemia": ["leuk"],
        "seizures": ["sz"],
        "sanfilippo": ["sanfilipo"],
        "sclerosis": ["schlerosis", "scleosis"],
        "respiratory": ["rep", "resp"],
        "syndrome": ["synd", "syndromes"],
        "undescended": ["undesc"],
    }
    replace_words = invert_dict_of_lists(REPLACE_ALTERNATES)

    for word, offset_start, _ in tokens_and_spans(cleaned_text, tokens):
        new_offset_start = offset_start
        unsplit_digits_and_spans = []  # Don't split up dates

        word = PRE_SPLIT_REPLACE.get(word, word)
        for w in word.split('/'):
            size = len(w)  # calculate on old size re-replace

            replace = replace_words.get(w.lower())
            if not replace and word.endswith("."):
                replace = replace_words.get(w[:-1].lower())

            if replace:
                w = replace

            new_offset_end = new_offset_start + size
            if w:
                data = (w, new_offset_start, new_offset_end)
                if w.isdigit():
                    unsplit_digits_and_spans.append(data)
                else:
                    handle_unsplit(unsplit_digits_and_spans, split_words_and_spans)
                    split_words_and_spans.append(data)
            new_offset_start = new_offset_end
        handle_unsplit(unsplit_digits_and_spans, split_words_and_spans)

    return split_words_and_spans


def sentences_and_offsets(txt):
    """ Split by nltk.sent_tokenize and new line """
    tokens = nltk.sent_tokenize(txt)

    split_sentences_and_offsets = []
    for sentence, sentence_offset, _ in tokens_and_spans(txt, tokens):
        new_offset = sentence_offset
        for line in sentence.split('\n'):
            if line:
                split_sentences_and_offsets.append((line, new_offset))
            new_offset += len(line) + 1

    return split_sentences_and_offsets


def process_text_phenotype(text_phenotype, hpo_records, hpo_word_lookup, hpo_single_words_by_length, omim_records, omim_word_lookup, omim_single_words_by_length, gene_records):
    tokenized_text_and_spans = word_tokenise_and_spans(text_phenotype.text)

    tokenized_text = []
    tokenized_spans = []
    for ts in tokenized_text_and_spans:
        tokenized_text.append(ts[0])
        tokenized_spans.append(ts[1:])

    tagged = nltk.pos_tag(tokenized_text)

    words = []
    tags = []
    spans = []
    commented_out = False
    for (w, t), span in zip(tagged, tokenized_spans):
        if not commented_out:
            if w == t:
                if w.startswith('--'):
                    commented_out = True
            else:
                words.append(w)
                tags.append(t)
                offset_start = span[0]
                offset_end = span[1]
                spans.append((offset_start, offset_end))

    words_and_spans = transform_words(words, tags, spans)

    try:
        parse_words(text_phenotype, words_and_spans, hpo_records, hpo_word_lookup, hpo_single_words_by_length, omim_records, omim_word_lookup, omim_single_words_by_length, gene_records)
    except SkipAllPhenotypeMatchException:
        logging.info("Completely skipping: %s", text_phenotype.text)

    text_phenotype.processed = True
    text_phenotype.save()


def break_up_hpo_terms(hpo_pks):
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


def get_omim_pks_by_term():
    """ Create entries for UNIQUE ';' separated terms
        ie "BECKWITH-WIEDEMANN SYNDROME; BWS" => "BECKWITH-WIEDEMANN" and "BECKWITH-WIEDEMANN SYNDROME" and "BWS" """

    omim_qs = OntologyTerm.objects.filter(ontology_service=OntologyService.OMIM)
    omim_pks_by_term = {}

    def break_up_dashes(omim_description, pk):
        if '-' in omim_description:
            cleaned_omim_description = omim_description.replace('-', ' ')
            omim_pks_by_term[cleaned_omim_description] = pk

    def break_up_syndromes_and_disease(omim_description, omim_alias):
        words = omim_description.split()
        if len(words) >= 2:
            for t in ['syndrome', 'disease']:
                try:
                    i = words.index(t)
                    if i >= 2:
                        term_before_syndrome = ' '.join(words[:i])
                        omim_pks_by_term[term_before_syndrome] = pk
                        break_up_dashes(term_before_syndrome, omim_alias)
                except ValueError:
                    pass

    for pk, name, aliases in omim_qs.values_list("pk", "name", "aliases"):
        for term in [name] + aliases:
            # Remove commas, as phenotype to match will have that done also
            term = term.strip().lower().replace(",", "")
            omim_pks_by_term[term] = pk
            for split_term in term.split(";"):
                omim_pks_by_term[split_term] = pk
                break_up_syndromes_and_disease(split_term, pk)
                break_up_dashes(split_term, pk)
    return omim_pks_by_term


def get_single_words_by_length(records, min_length):
    """ if key is pure alpha (ie FOO not FOO,BAR or FOO1 or FOO BAR) put lookup into dict by lengths """
    words_by_length = defaultdict(dict)
    for k, v in records.items():
        if k.isalpha():
            length = len(k)
            if length >= min_length:
                words_by_length[length][k] = v
    return words_by_length


def default_lookup_factory():
    # Start with synonyms so they're overwritten by HPO
    hpo_qs = OntologyTerm.objects.filter(ontology_service=OntologyService.HPO)
    hpo_pks = {}
    for pk, name, aliases in hpo_qs.values_list('pk', 'name', "aliases"):
        for alias in [name] + aliases:
            k = alias.lower().replace(",", "")
            hpo_pks[k] = pk

    break_up_hpo_terms(hpo_pks)
    hpo_word_lookup = create_word_lookups(hpo_pks)

    omim_pks = get_omim_pks_by_term()
    omim_word_lookup = create_word_lookups(omim_pks)

    gene_symbol_records = dict(GeneSymbol.objects.annotate(lower=Lower("pk")).values_list("lower", "pk"))
    hpo_single_words_by_length = get_single_words_by_length(hpo_pks, 5)
    omim_single_words_by_length = get_single_words_by_length(omim_pks, 5)
    return hpo_pks, hpo_word_lookup, hpo_single_words_by_length, omim_pks, omim_word_lookup, omim_single_words_by_length, gene_symbol_records


def cached_lookup_factory(*args):
    """ Use like so:
            lookups = default_lookup_factory() # Keep in scope
            lookup_factory = cached_lookup_factory(*lookups)
    """

    def lookup_factory():
        return args

    return lookup_factory


def replace_comments_with_spaces(text):
    """ Replace with spaces instead of removing so that characters offsets remain the same """
    COMMENT = '-'
    COMMENT_REPLACE = ' '
    cleaned_chars = []

    last_char = None
    in_comment = False
    for c in text:
        if c == COMMENT and last_char == COMMENT:
            in_comment = True
        elif c == '\n':
            in_comment = False

        if in_comment:
            cleaned_chars.append(COMMENT_REPLACE)
        else:
            cleaned_chars.append(c)
        last_char = c

    return ''.join(cleaned_chars)


def create_phenotype_description(text, lookup_factory=default_lookup_factory):
    phenotype_description = PhenotypeDescription.objects.create(original_text=text)
    known_sentences = []
    unknown_sentences = []

    # Strip carriage returns - some browsers send \r\n in forms, Ajax usually sends \n
    # These need to be the same so save/ajax preview look the same
    text = text.replace('\r', '')
    text = replace_comments_with_spaces(text)  # Need to keep length the same for offsets
    for sentence, sentence_offset in sentences_and_offsets(text):
        text_phenotype, tp_created = TextPhenotype.objects.get_or_create(text=sentence)
        tps = TextPhenotypeSentence.objects.create(phenotype_description=phenotype_description,
                                                   text_phenotype=text_phenotype,
                                                   sentence_offset=sentence_offset)
        if tp_created:
            has_alpha_numeric = any([x.isalnum() for x in sentence])
            if not has_alpha_numeric:  # Impossible to match anything so skip
                text_phenotype.processed = True
                text_phenotype.save()
                known_sentences.append(tps)
            else:
                unknown_sentences.append(tps)
        else:
            known_sentences.append(tps)

    if unknown_sentences:
        args = lookup_factory()
        for sentence in unknown_sentences:
            process_text_phenotype(sentence.text_phenotype, *args)
            known_sentences.append(sentence)

    return phenotype_description


def bulk_patient_phenotype_matching(patients=None):
    if patients is None:
        patients = Patient.objects.filter(phenotype__isnull=False)

    start = time.time()
    lookups = default_lookup_factory()
    get_and_log_time_since(start, "load references")

    lookup_factory = cached_lookup_factory(*lookups)

    start = time.time()
    patients = list(patients)
    num_patients = len(patients)
    num_parsed_phenotypes = 0
    if num_patients:
        for i, p in enumerate(patients):

            parsed_phenotypes = p.process_phenotype_if_changed(lookup_factory=lookup_factory)
            num_parsed_phenotypes += parsed_phenotypes
            if not i % 50:
                perc_complete = 100.0 * i / num_patients
                logging.info("%d patients %0.2f%% complete", i, perc_complete)

        ts = get_and_log_time_since(start, "bulk_patient_phenotype_matching")
        if num_parsed_phenotypes:
            time_per_patient = ts / num_parsed_phenotypes
            pps = 1.0 / time_per_patient
            logging.info("%d parsed patient phenotypes - %.2f parsed per second/%.2f seconds per patient", num_parsed_phenotypes, pps, time_per_patient)

            total_matches = TextPhenotypeMatch.objects.count()
            tpm_qs = TextPhenotypeMatch.objects.values("match_type")
            tpm_qs = tpm_qs.annotate(count=Count("pk")).values_list("match_type", "count")
            logging.info("%d match types:", total_matches)
            for match_type, count in tpm_qs:
                logging.info("%s: %d", match_type, count)
    else:
        logging.info("No patients")
