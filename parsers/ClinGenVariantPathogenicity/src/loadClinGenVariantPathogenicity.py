import os
import csv
from Common.extractor import Extractor
from Common.loader_interface import SourceDataLoader
from biolink_constants import PRIMARY_KNOWLEDGE_SOURCE, NODE_TYPES, SEQUENCE_VARIANT
from Common.prefixes import HGNC  # only an example, use existing curie prefixes or add your own to the prefixes file
from Common.utils import GetData
from Common.utils import LoggingUtil
from datetime import date


##############
# Class: ClinGenVariantPathogenicity  source loader
# Desc: Class that loads/parses the ClinGenVariantPathogenicity data.
##############
class ClinGenVariantPathogenicityLoader(SourceDataLoader):
    source_id: str = "ClinGenVariantPathogenicity"
    provenance_id: str = "infores:clingen"
    # increment parsing_version whenever changes are made to the parser that would result in changes to parsing output
    parsing_version: str = "1.0"
    has_sequence_variants = (
        True  # Flag to use robokop_genetics server to tackle sequence variant data
    )

    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode - sets the run into test mode
        :param source_data_dir - the specific storage directory to save files in
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)
        self.data_url = "http://erepo.clinicalgenome.org/evrepo/api/classifications/"
        self.clingen_variant_pathogenicity_file = "all?format=tabbed"
        self.data_files = [self.clingen_variant_pathogenicity_file]

    def get_latest_source_version(self) -> str:
        # No version is available at the source, using the year_month when the code was run as versioning proxy
        latest_version = date.today().strftime("%Y%m")
        return latest_version

    def get_data(self) -> bool:
        # get_data is responsible for fetching the files in self.data_files and saving them to self.data_path
        source_data_url = f"{self.data_url}{self.clingen_variant_pathogenicity_file}"
        data_puller = GetData()
        data_puller.pull_via_http(source_data_url, self.data_path)
        return True

    def parse_data(self) -> dict:
        """
        Parses the data file for graph nodes/edges

        :return: ret_val: load_metadata
        """
        extractor = Extractor(file_writer=self.output_file_writer)
        clingen_variant_pathogenicity_file: str = os.path.join(
            self.data_path, self.clingen_variant_pathogenicity_file
        )

        # TODO: encourage upstream to pass extra args to parse_row so we can exclude_unconnected_nodes
        with open(clingen_variant_pathogenicity_file, "rt") as fp:
            reader = csv.DictReader(fp, dialect="excel-tab")
            extractor.json_extract(
                reader,
                lambda line: f"CLINVARVARIANT:{line['ClinVar Variation Id']}",  # subject id
                lambda line: line["Mondo Id"],  # object id
                lambda line: (
                    "causes" if line["Retracted"] == "false" else None
                ),  # predicate extractor
                lambda line: {
                    NODE_TYPES: SEQUENCE_VARIANT,
                    "Variation": line["#Variation"],
                    "HGNC_Gene_Symbol": line["HGNC Gene Symbol"],
                },  # subject properties
                lambda line: {},  # object properties
                lambda line: {
                    PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                    "Assertion": line["Assertion"],
                    "Mode_Of_Inheritance": moi_normalizer(
                        line["Mode of Inheritance"],
                        line["Evidence Repo Link"],
                    ),
                    "Applied_Evidence_Codes_Met": line["Applied Evidence Codes (Met)"],
                    "Applied_Evidence_Codes_Not_Met": line[
                        "Applied Evidence Codes (Not Met)"
                    ],
                    "Summary": line["Summary of interpretation"],
                    "Pubmed_Articles": line["PubMed Articles"],
                    "Expert_Panel": line["Expert Panel"],
                    "Evidence_Repo_Link": line["Evidence Repo Link"],
                    "Guideline": line["Guideline"],
                    "Approval_Date": line["Approval Date"],
                    "Published_Date": line["Published Date"],
                    **get_edge_properties(line["Assertion"]),
                },  # edge properties
            )
        return extractor.load_metadata


# Supporting functions for processing the data

# Logging
logger = LoggingUtil.init_logging(
    "ORION.parsers.ClinGenVariantPathogenicityLoader",
    log_file_path=os.environ["ORION_LOGS"],
)

moi_lookup = {
    "Autosomal dominant inheritance": "HP:0000006",
    "Autosomal dominant inheritance (with paternal imprinting (HP:0012274))": "HP:0012274",
    "Autosomal dominant inheritance (mosaic)": ["HP:0000006", "HP:0001442"],
    "Autosomal recessive inheritance": "HP:0000007",
    "Autosomal recessive inheritance (with genetic anticipation)": ["HP:0000007"],
    "X-linked inheritance": "HP:0001417",
    "X-linked inheritance (dominant (HP:0001423))": "HP:0001423",
    "X-linked inheritance (recessive (HP:0001419))": "HP:0001419",
    "Semidominant inheritance": "HP:0032113",
    "Mitochondrial inheritance": "HP:0001427",
    "Mitochondrial inheritance (primarily or exclusively heteroplasmic)": "HP:0001427",  # No HPO term for heteroplasmic type
}


def moi_normalizer(MOI, EREPO_LINK):
    MOI = str(MOI)
    try:
        HPO = moi_lookup[MOI]
    except KeyError:
        logger.warning(
            f"We do not have a mapping for {MOI=} in the moi_lookup dictionary at source {EREPO_LINK}"
        )
        HPO = ""
    # TODO: Check the second HPO for "Autosomal recessive inheritance (with genetic anticipation)" and "Mitochondrial inheritance (primarily or exclusively heteroplasmic)"
    return {"label": MOI, "HPO": HPO}


def get_edge_properties(assertion):
    if assertion == "Benign" or assertion == "Likely Benign":
        return {"direction": "Contradicts", "negated": True}
    elif assertion == "Likely Pathogenic" or assertion == "Pathogenic":
        return {"direction": "Supports", "negated": False}
    elif assertion == "Uncertain Significance":
        return {"direction": "Inconclusive", "negated": True}
    else:
        return {"Status": "Not evaluated", "direction": "Inconclusive", "negated": True}
