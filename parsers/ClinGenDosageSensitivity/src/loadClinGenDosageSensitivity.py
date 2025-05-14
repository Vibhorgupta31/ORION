import os
import enum
import gzip
from Common.extractor import Extractor
from Common.loader_interface import SourceDataLoader
from biolink_constants import PRIMARY_KNOWLEDGE_SOURCE, NODE_TYPES, SEQUENCE_VARIANT
from Common.utils import GetData
from datetime import date


# For parsing tsv files
# HI: HaploInsufficiency ; TS : TriploSensitivity
class ClinGenDosageSensitivityCOLS(enum.IntEnum):
    REGION = 0
    GENE = 1
    HI_DISEASE = -2
    TS_DISEASE = -1
    HI_SCORE = 4
    HI_DESCRIPTION = 5
    TS_SCORE = 12
    TS_DESCRIPTION = 13


HUMAN_DISEASE = "MONDO:0700096"


##############
# Class: ClinGenDosageSensitivity  data loader
# Desc: Class that loads/parses the ClinGen Dosage Sensitivity  data.
###########
class ClinGenDosageSensitivityLoader(SourceDataLoader):
    source_id: str = "ClinGenDosageSensitivity"
    provenance_id: str = "infores:clingen"
    # increment parsing_version whenever changes are made to the parser that would result in changes to parsing output
    parsing_version: str = "v1.0"

    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode - sets the run into test mode
        :param source_data_dir - the specific storage directory to save files in
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)
        self.cligen_dosage_sensitivity_url = "http://ftp.clinicalgenome.org/"
        self.cligen_dosage_sensitivity_gene_file = (
            "ClinGen_gene_curation_list_GRCh38.tsv"
        )
        self.clingen_dosage_sensitivity_region_file = (
            "ClinGen_region_curation_list_GRCh38.tsv"
        )
        self.data_files = [
            self.cligen_dosage_sensitivity_gene_file,
            self.clingen_dosage_sensitivity_region_file,
        ]

    def get_latest_source_version(self) -> str:
        # No version is available at the source, using the year_month when the code was run as versioning proxy
        latest_version = date.today().strftime("%Y%m")
        return latest_version

    def get_data(self) -> bool:
        # get_data is responsible for fetching the files in self.data_files and saving them to self.data_path
        data_puller = GetData()
        for source in self.data_files:
            source_data_url = f"{self.cligen_dosage_sensitivity_url}{source}"
            data_puller.pull_via_http(source_data_url, self.data_path)
        return True

    def dosage_sensitivity_edge_generator(
        self, data_file: str, subject_extractor, predicate
    ):
        """
        Generator function to yield edges from the dosage sensitivity data file.

        :param data_file: Path to the data file.
        :param subject_extractor: Function to extract subject ID.
        :param predicate: predicate associated with the relationship between subject and object.
        :return: Generator yielding edges as dictionaries.
        """
        with open(data_file, "rt") as fp:
            for line in fp:
                if line.startswith("#"):
                    continue
                if len(line) == 0:
                    continue
                line = line.strip("\n").split("\t")

                record = {
                    "subject": subject_extractor(line),
                    "object": line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value]
                    or HUMAN_DISEASE,
                    "predicate": predicate,
                    "subject_properties": {},
                    "object_properties": {},
                    "edge_properties": {
                        PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                        "Haploinsufficiency Description": line[
                            ClinGenDosageSensitivityCOLS.HI_DESCRIPTION.value
                        ],
                        "Haploinsufficiency Score": line[
                            ClinGenDosageSensitivityCOLS.HI_SCORE.value
                        ],
                        **get_edge_properties(
                            line[ClinGenDosageSensitivityCOLS.HI_SCORE.value],
                            line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value],
                        ),
                    },
                }
                if (
                    line[ClinGenDosageSensitivityCOLS.HI_SCORE.value]!= "Not yet evaluated"):
                    yield record

                record["object"] = (
                    line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value] or HUMAN_DISEASE
                )
                record["edge_properties"] = {
                    PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                    "Triplosensitivity Description": line[
                        ClinGenDosageSensitivityCOLS.TS_DESCRIPTION.value
                    ],
                    "Triplosensitivity Score": line[
                        ClinGenDosageSensitivityCOLS.TS_SCORE.value
                    ],
                    **get_edge_properties(
                        line[ClinGenDosageSensitivityCOLS.TS_SCORE.value],
                        line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value],
                    ),
                }
                if (line[ClinGenDosageSensitivityCOLS.TS_SCORE.value]!= "Not yet evaluated"):
                    yield record

    def parse_data(self) -> dict:
        """
        Parses the data file for graph nodes/edges

        :return: ret_val: load_metadata
        """
        extractor = Extractor(file_writer=self.output_file_writer)
        # Gene file processing
        dosage_sensitivity_gene_file: str = os.path.join(
            self.data_path, self.cligen_dosage_sensitivity_gene_file
        )
        extractor.json_extract(
            self.dosage_sensitivity_edge_generator(
                dosage_sensitivity_gene_file,
                subject_extractor=lambda line: "NCBIGene:%s"
                % line[ClinGenDosageSensitivityCOLS.GENE.value],
                predicate="gene associated with condition",
            )
        )
        # Region file processing
        dosage_sensitivity_region_file: str = os.path.join(
            self.data_path, self.clingen_dosage_sensitivity_region_file
        )
        extractor.json_extract(
            self.dosage_sensitivity_edge_generator(
                dosage_sensitivity_region_file,
                subject_extractor=lambda line: line[
                    ClinGenDosageSensitivityCOLS.REGION.value
                ],
                predicate="region associated with condition",
            )
        )

        return extractor.load_metadata


# Created a common function to take the score value and return attributes,
# this may have problems with TS and HI score where the mode of inheritance is captured, followed ClinGen criteria for converting scores
# Numeric scores are as describe in the
# [ClinGen Dosage Sensitivity Curation Guidelines](https://clinicalgenome.org/docs/dosage-standard-operating-procedure-scoring-guide)
# 0: No evidence available
# 1: Little evidence for dosage pathogenicity
# 2: Some evidence for dosage pathogenicity
# 3: Sufficient evidence for dosage pathogenicity
# 30: Gene associated with autosomal recessive phenotype (so haploinsufficiency not applicable)
# 40: Dosage sensitivity unlikely
def get_edge_properties(score, disease_id):
    if disease_id != "":
        try:
            score = int(score)
        except ValueError:
            return {"Status": "Not yet evaluated"}
        if score in (2, 3):
            return {"negated": False}
        elif score == 30:
            return {"negated": False}
        elif score in (0, 1, 40):  #  condition for negation
            return {"negated": True}
        elif score == -1:
            return {"Status": "Not planned to be evaluated", "negated": True}
    else:  # Negate only if no HI disease and TS disease is available for the gene
        return {"negated": True}
