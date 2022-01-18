from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Tuple
from classification.models import ClassificationModification
from classification.views.exports.classification_export_filter import AlleleData, ClassificationFilter
from genes.hgvs import CHGVS


@dataclass
class CHGVSData:
    """
    A sub-division of AlleleData.
    Will create one record per unique c.hgvs string within the allele*
    (c.hgvs differing in just transcript version are still bundled together)

    :var allele: The allele data record
    :var chgvs: The c.hgvs with the highest found transcript version
    :var different_chgvs: Bool indicating if multiple c.hgvs versions were bundled together here
    :var cms: The classifications
    """
    allele: AlleleData

    @property
    def source(self) -> ClassificationFilter:
        return self.allele.source

    chgvs: CHGVS
    different_chgvs: bool = False
    cms: List[ClassificationModification] = field(default_factory=list)

    @staticmethod
    def split_into_c_hgvs(
            allele_data: AlleleData,
            use_full: bool) -> List['CHGVSData']:
        """
        Break up an AlleleData into sub CHGVSDatas
        :param allele_data: The Alissa data to split
        :param use_full: Should the c.hgvs use explicit bases when optional (required by Alissa)
        :return: An array of c.hgvs based data, most often will only be 1 record
        """
        genome_build = allele_data.source.genome_build
        cms = allele_data.cms
        by_chgvs: Dict[CHGVS, List[Tuple[ClassificationModification, int]]] = defaultdict(list)
        for cm in cms:
            c_hgvs = CHGVS(cm.classification.get_c_hgvs(genome_build=genome_build, use_full=use_full))
            by_chgvs[c_hgvs.without_transcript_version].append((cm, c_hgvs.transcript_parts.version))

        sub_datas: List[CHGVSData] = list()
        for c_hgvs, versioned_cms in by_chgvs.items():
            cms = [vcm[0] for vcm in versioned_cms]
            versions = [v for v in reversed(sorted(set(vcm[1] for vcm in versioned_cms))) if v is not None]
            # always provide highest transcript version when there's multiple
            # but make sure we at least have 1 transcript version and not an LRG for example
            if versions:
                c_hgvs.transcript = c_hgvs.transcript + f'.{versions[0]}'

            sub_datas.append(CHGVSData(allele=allele_data, chgvs=c_hgvs, different_chgvs=False, cms=cms))

        # If we actually want to merge c.hgvs data together
        if len(sub_datas) > 1:
            first = sorted(sub_datas, key=lambda cd: (cd.chgvs.transcript_parts.version, cd.chgvs.sort_str), reverse=True)[0]
            first.cms = allele_data.cms
            first.different_chgvs = True
            return [first]

        return sub_datas
