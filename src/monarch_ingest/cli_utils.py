import csv
import os
import sys
import tarfile

import duckdb
import yaml
from pathlib import Path
from typing import Optional
from biolink_model.datamodel import model  # import the pythongen biolink model to get the version
from linkml_runtime import SchemaView
from linkml.utils.helpers import convert_to_snake_case
import requests

# from loguru import logger
import pandas
import sh

from cat_merge.merge import merge
from closurizer.closurizer import add_closure
from kgx.cli.cli_utils import transform as kgx_transform
from koza.cli_utils import transform_source
from koza.model.config.source_config import OutputFormat
from linkml_runtime.utils.formatutils import camelcase

from monarch_ingest.utils.ingest_utils import ingest_output_exists, file_exists, get_ingests
from monarch_ingest.utils.log_utils import get_logger
from monarch_ingest.utils.export_utils import export


OUTPUT_DIR = "output"


def transform_one(
    ingest: Optional[str] = None,
    output_dir: str = OUTPUT_DIR,
    row_limit: Optional[int] = None,
    rdf: bool = False,
    force: Optional[bool] = False,
    verbose: Optional[bool] = None,
    log: bool = False,
):
    logger = get_logger(name=ingest if log else None, verbose=verbose)

    ingests = get_ingests()

    if ingest not in ingests:
        # if log: logger.removeHandler(fh)
        raise ValueError(f"{ingest} is not a valid ingest - see ingests.yaml for a list of options")

    # if a url is provided instead of a config, just download the file and copy it to the output dir
    if "url" in ingests[ingest]:
        for url in ingests[ingest]["url"]:
            filename = url.split("/")[-1]

            if Path(f"{output_dir}/transform_output/{filename}").is_file() and not force:
                continue

            response = requests.get(url, allow_redirects=True)
            with open(f"{output_dir}/transform_output/{filename}", "wb") as f:
                f.write(response.content)
        return

    source_file = Path(Path(__file__).parent, ingests[ingest]["config"])

    if not Path(source_file).is_file():
        # if log: logger.removeHandler(fh)
        raise ValueError(f"Source file {source_file} does not exist")

    if ingest_output_exists(ingest, output_dir) and not force:
        logger.info(f"Transformed output exists - skipping ingest: {ingest} - To run this ingest anyway, use --force")
        # if log: logger.removeHandler(fh)
        return

    logger.info(f"Running ingest: {ingest}")
    try:
        transform_source(
            source=source_file.as_posix(),
            output_dir=f"{output_dir}/transform_output",
            output_format=OutputFormat.tsv,
            row_limit=row_limit,
            verbose=verbose,
            # verbose=False if verbose == None else verbose,
            # log=log
        )
    except ValueError as e:
        # if log: logger.removeHandler(fh)
        raise ValueError(f"Missing data - {e}")

    if rdf:
        logger.info(f"Creating rdf output {output_dir}/rdf/{ingest}.nt.gz ...")

        Path(f"{output_dir}/rdf").mkdir(parents=True, exist_ok=True)

        src_files = []
        src_nodes = f"{output_dir}/transform_output/{ingest}_nodes.tsv"
        src_edges = f"{output_dir}/transform_output/{ingest}_edges.tsv"
        if Path(src_nodes).is_file():
            src_files.append(src_nodes)
        if Path(src_edges).is_file():
            src_files.append(src_edges)

        kgx_transform(
            inputs=src_files,
            input_format="tsv",
            stream=True,
            output=f"{output_dir}/rdf/{ingest}.nt.gz",
            output_format="nt",
            output_compression="gz",
        )

    if not ingest_output_exists(ingest, f"{output_dir}/transform_output"):
        # if log: logger.removeHandler(fh)
        raise FileNotFoundError(f"Ingest {ingest} did not produce the the expected output")

    # if log: logger.removeHandler(fh)


def transform_phenio(
    output_dir: str = OUTPUT_DIR,
    force: bool = False,
    verbose: Optional[bool] = False,
    log: bool = False,
):
    # if log: fh = add_log_fh(logger, "logs/phenio.log")
    logger = get_logger(name="phenio" if log else None, verbose=verbose)

    # TODO: can this be fetched from a purl?
    biolink_model_schema = SchemaView(
        f"https://raw.githubusercontent.com/biolink/biolink-model/v{model.version}/biolink-model.yaml"
    )

    phenio_tar = "data/monarch/kg-phenio.tar.gz"
    if not Path(phenio_tar).is_file():
        # if log: logger.removeHandler(fh)
        raise FileNotFoundError("data/monarch/kg-phenio.tar.gz")

    nodefile = "merged-kg_nodes.tsv"
    edgefile = "merged-kg_edges.tsv"

    tar = tarfile.open(phenio_tar)
    tar.extract(nodefile, "data/monarch")
    tar.extract(edgefile, "data/monarch")

    Path(f"{output_dir}/transform_output").mkdir(parents=True, exist_ok=True)

    nodes = f"{output_dir}/transform_output/phenio_nodes.tsv"
    edges = f"{output_dir}/transform_output/phenio_edges.tsv"

    if (force is False) and file_exists(nodes) and file_exists(edges):
        logger.info(f"Transformed output exists - skipping ingest: Phenio - To run this ingest anyway, use --force")
        # if log: logger.removeHandler(fh)
        return

    nodes_df = pandas.read_csv(
        f"data/monarch/{nodefile}", sep="\t", dtype="string", quoting=csv.QUOTE_NONE, lineterminator="\n"
    )
    nodes_df = nodes_df[~nodes_df["id"].str.contains("omim.org|hgnc_id")]
    nodes_df = nodes_df[~nodes_df["id"].str.startswith("MGI:")]

    # Hopefully this won't be necessary long term, but these IDs are coming
    # in with odd OBO prefixes from Phenio currently.
    obo_prefixes_to_repair = ["FBbt", "WBbt", "ZFA", "XAO"]
    for prefix in obo_prefixes_to_repair:
        nodes_df["id"] = nodes_df["id"].str.replace(f"OBO:{prefix}_", f"{prefix}:")

    # These bring in nodes necessary for other ingests, but won't capture the same_as / equivalentClass
    # associations that we'll also need
    exclude_prefixes = ["HGNC", "FlyBase", "http", "biolink"]

    Path(f"{output_dir}/qc/").mkdir(parents=True, exist_ok=True)
    excluded_nodes = nodes_df[nodes_df["id"].str.startswith(tuple(exclude_prefixes))]
    nodes_df = nodes_df[~nodes_df["id"].str.startswith(tuple(exclude_prefixes))]

    # Replace biolink:Occurrent category with biolink:BiologicalProcessOrActivity as a fallback to rescue
    # GO nodes that are getting a mixin category of Occurrent that we don't want to exclude all of
    nodes_df["category"] = nodes_df["category"].str.replace("biolink:Occurrent", "biolink:BiologicalProcessOrActivity")

    valid_node_categories = {
        f"biolink:{camelcase(cat)}" for cat in biolink_model_schema.class_descendants("named thing")
    }
    phenio_node_categories = set(nodes_df["category"].unique())
    invalid_node_categories = phenio_node_categories - valid_node_categories
    if invalid_node_categories:
        logger.error(f"Invalid node categories: {invalid_node_categories}")
        invalid_node_categories_df = nodes_df[nodes_df["category"].isin(invalid_node_categories)]
        excluded_nodes = pandas.concat([excluded_nodes, invalid_node_categories_df])
        logger.error(f"Removing {len(invalid_node_categories_df)} nodes with invalid categories")
        nodes_df = nodes_df[~nodes_df["category"].isin(invalid_node_categories)]

    excluded_nodes.to_csv(f"{output_dir}/qc/excluded_phenio_nodes.tsv", sep="\t", index=False)

    nodes_df.to_csv(nodes, sep="\t", index=False)

    edges_df = pandas.read_csv(
        f"data/monarch/{edgefile}", sep="\t", dtype="string", quoting=csv.QUOTE_NONE, lineterminator="\n"
    )
    edges_df.drop(
        edges_df.columns.difference(
            [
                "id",
                "subject",
                "predicate",
                "object",
                "category",
                "primary_knowledge_source",
                "aggregator_knowledge_source",
                "knowledge_level",
                "agent_type",
            ]
        ),
        axis=1,
        inplace=True,
    )

    # if knowledge level doesn't exist, add it and assign to knowledge_assertion
    if "knowledge_level" not in edges_df.columns:
        edges_df["knowledge_level"] = "knowledge_assertion"
    # same for agent_type, setting it to manual_agent
    if "agent_type" not in edges_df.columns:
        edges_df["agent_type"] = "manual_agent"

    # prepend infores:monarchinitiative to the aggregator_knowledge_source column for edges that don't have it
    edges_df["aggregator_knowledge_source"] = edges_df["aggregator_knowledge_source"].apply(
        lambda x: f"infores:monarchinitiative|{x}" if not x.startswith("infores:monarchinitiative") else x
    )

    edges_df = edges_df[edges_df["predicate"].str.contains(":")]

    # assign level association category if edge category is empty
    edges_df["category"].fillna("biolink:Association", inplace=True)

    # Hopefully this won't be necessary long term, but these IDs are coming
    # in with odd OBO prefixes from Phenio currently.
    for prefix in obo_prefixes_to_repair:
        for field in ["subject", "object"]:
            edges_df[field] = edges_df[field].str.replace(f"OBO:{prefix}_", f"{prefix}:")

    # Only keep edges where the subject and object both are within our allowable prefix list
    edges_df = edges_df[
        ~edges_df["subject"].str.startswith(tuple(exclude_prefixes))
        & ~edges_df["object"].str.startswith(tuple(exclude_prefixes))
    ]

    # Remove edges where the subject or object is NA
    edges_df = edges_df[~edges_df["subject"].isna() & ~edges_df["object"].isna()]

    valid_predicates = {
        f"biolink:{convert_to_snake_case(pred)}" for pred in biolink_model_schema.slot_descendants("related to")
    }
    phenio_predicates = set(edges_df["predicate"].unique())
    invalid_predicates = phenio_predicates - valid_predicates
    if invalid_predicates:
        logger.error(f"Invalid predicates found in Phenio associations: {invalid_predicates}")
        invalid_edges_df = edges_df[edges_df["predicate"].isin(invalid_predicates)]
        logger.error(f"Removing {invalid_edges_df.shape[0]} edges with invalid predicates")
        edges_df = edges_df[~edges_df["predicate"].isin(invalid_predicates)]

    valid_edge_categories = {
        f"biolink:{camelcase(cat)}" for cat in biolink_model_schema.class_descendants("association")
    }
    phenio_edge_categories = set(edges_df["category"].unique())
    invalid_edge_categories = phenio_edge_categories - valid_edge_categories
    if invalid_edge_categories:
        logger.error(f"Invalid edge categories: {invalid_edge_categories}")
        invalid_edge_categories_df = edges_df[edges_df["category"].isin(invalid_edge_categories)]
        logger.error(f"Removing {len(invalid_edge_categories_df)} edges with invalid categories")
        edges_df = edges_df[~edges_df["category"].isin(invalid_edge_categories)]

    edges_df.to_csv(edges, sep="\t", index=False)
    Path(f"data/monarch/{nodefile}").unlink()
    Path(f"data/monarch/{edgefile}").unlink()

    if not file_exists(nodes) or not file_exists(edges):
        # if log: logger.removeHandler(fh)
        raise FileNotFoundError("Phenio transform did not produce the expected output")

    # if log: logger.removeHandler(fh)


def transform_all(
    output_dir: str = OUTPUT_DIR,
    row_limit: Optional[int] = None,
    rdf: bool = False,
    force: bool = False,
    verbose: Optional[bool] = None,
    log: bool = False,
):
    # if log: fh = add_log_fh(logger, Path(f"logs/all_ingests.log"))
    logger = get_logger(name="all_ingests" if log else None, verbose=verbose)

    # TODO:
    # - check for data - download if missing (maybe y/n prompt?)
    # - check for difference in data? maybe implement in kghub downloder instead?

    try:
        transform_phenio(output_dir=output_dir, force=force)
    except Exception as e:
        logger.error(f"Error running Phenio ingest: {e}")

    ingests = get_ingests()
    for ingest in ingests:
        try:
            transform_one(
                ingest=ingest,
                output_dir=output_dir,
                row_limit=row_limit,
                rdf=rdf,
                force=force,
                verbose=verbose,
                log=log,
            )
        except Exception as e:
            logger.error(f"Error running ingest {ingest}: {e}")
            pass

    # if log: logger.removeHandler(fh)


def get_release_version():
    import datetime

    return datetime.datetime.now().strftime("%Y-%m-%d")


def get_data_versions(output_dir: str = OUTPUT_DIR):
    import requests as r

    data = {}
    data["phenio"] = r.get("https://api.github.com/repos/monarch-initiative/phenio/releases").json()[0]["tag_name"]
    data["alliance"] = r.get("https://fms.alliancegenome.org/api/releaseversion/current").json()["releaseVersion"]
    Path(f"{output_dir}").mkdir(parents=True, exist_ok=True)
    with open(f"{output_dir}/metadata.yaml", "w") as f:
        f.write("data:\n")
        for data_source, version in data.items():
            f.write(f"  {data_source}: {str(version)}\n")


def get_pkg_versions(output_dir: str = OUTPUT_DIR, release_version: Optional[str] = None):
    from importlib.metadata import version

    packages = {}
    packages["biolink"] = version("biolink-model")
    packages["koza"] = version("koza")
    packages["monarch-ingest"] = version("monarch-ingest")
    kg_version = get_release_version() if release_version is None else release_version
    with open("data/metadata.yaml", "r") as f:
        data_versions = yaml.load(f, Loader=yaml.FullLoader)["data"]

    Path(f"{output_dir}").mkdir(parents=True, exist_ok=True)
    with open(f"{output_dir}/metadata.yaml", "w") as f:
        f.write(f"kg-version: {str(kg_version)}\n\n")
        f.write("packages:\n")
        for k, v in packages.items():
            f.write(f"  {k}: {str(v)}\n")
        f.write("\ndata:\n")
        for k, v in data_versions.items():
            f.write(f"  {k}: {str(v)}\n")


def merge_files(
    name: str = "monarch-kg",
    input_dir: str = f"{OUTPUT_DIR}/transform_output",
    output_dir: str = OUTPUT_DIR,
    verbose: Optional[bool] = None,
):
    logger = get_logger(None, verbose)
    logger.info("Generating mappings...")

    mappings = []
    mappings.append("data/monarch/mondo.sssom.tsv")
    mappings.append("data/monarch/gene_mappings.sssom.tsv")
    mappings.append("data/monarch/mesh_chebi_biomappings.sssom.tsv")

    logger.info("Merging knowledge graph...")

    merge(name=name, source=input_dir, output_dir=output_dir, mappings=mappings)


def apply_closure(
    name: str = "monarch-kg",
    closure_file: str = f"data/monarch/phenio-relation-filtered.tsv",
    output_dir: str = OUTPUT_DIR,
):
    edges_output_file = f"{output_dir}/{name}-denormalized-edges.tsv"
    nodes_output_file = f"{output_dir}/{name}-denormalized-nodes.tsv"
    database = f"{name}.duckdb"
    add_closure(
        kg_archive=f"{output_dir}/{name}.tar.gz",
        closure_file=closure_file,
        edges_output_file=edges_output_file,
        nodes_output_file=nodes_output_file,
        # full closure expansion for these fields
        edge_fields=[
            "subject",
            "object",
            "disease_context_qualifier",
        ],
        #  just populate _label, _category, _namespace properties for these fields
        edge_fields_to_label=[
            "qualifiers",
            "species_context_qualifier",
            "stage_qualifier",
            "sex_qualifier",
            "onset_qualifier",
            "frequency_qualifier",
        ],
        node_fields=["has_phenotype"],
        evidence_fields=["has_evidence", "publications"],
        additional_node_constraints="has_phenotype_edges.negated is null or has_phenotype_edges.negated = 'False'",
        grouping_fields=["subject", "negated", "predicate", "object"],
    )
    sh.mv(database, f"{output_dir}/")


def load_sqlite():
    sh.bash("scripts/load_sqlite.sh")


def load_solr():
    sh.bash("scripts/load_solr.sh", _out=sys.stdout, _err=sys.stderr)


def load_jsonl():
    db = duckdb.connect('output/monarch-kg.duckdb')

    biolink_model = SchemaView(
        f"https://raw.githubusercontent.com/biolink/biolink-model/v{model.version}/biolink-model.yaml"
    )

    class_ancestor_dict = {}
    for c in biolink_model.all_classes().values():
        class_ancestor_dict[f"biolink:{camelcase(c.name)}"] = [
            f"biolink:{camelcase(a)}" for a in biolink_model.class_ancestors(c.name)
        ]
    all_slot_names = biolink_model.all_slots().keys()

    # this may appear to be unused, but it's accessed in the duckdb sql queries below
    class_ancestor_df = pandas.DataFrame(list(class_ancestor_dict.items()), columns=['classname', 'ancestors'])

    node_columns = db.sql("PRAGMA table_info(nodes);").df()["name"].to_list()
    edge_columns = db.sql("PRAGMA table_info(edges);").df()["name"].to_list()

    def slot_is_multi_valued(slot_name: str) -> bool:
        slot_name = slot_name.lower().replace("_", " ")
        if slot_name not in all_slot_names:
            return False
        return biolink_model.get_slot(slot_name).multivalued

    mv_node_columns = [col for col in node_columns if slot_is_multi_valued(col) and col != "category"]
    mv_edge_columns = [col for col in edge_columns if slot_is_multi_valued(col) and col != "category"]
    mv_node_replacement = ", ".join([f"string_split({col}, '|') as {col}" for col in mv_node_columns])
    mv_edge_replacement = ", ".join([f"string_split({col}, '|') as {col}" for col in mv_edge_columns])

    db.sql(
        f"""
    copy (
      select nodes.* replace (ancestors as category, {mv_node_replacement}) 
      from nodes
        join class_ancestor_df on category = classname  
    ) to 'output/monarch-kg_nodes.jsonl' (FORMAT JSON);
    """
    )

    db.sql(
        f"""
    copy (
      select edges.* replace (ancestors as category, {mv_edge_replacement}),  
      from edges
        join class_ancestor_df on category = classname  
    ) to 'output/monarch-kg_edges.jsonl' (FORMAT JSON);
    """
    )


def create_qc_reports():
    database_file = "output/monarch-kg.duckdb"
    # error if the database exists but needs to be gunzipped
    if Path(database_file + ".gz").is_file():
        raise FileExistsError(database_file + ".gz", "Database exists but needs to be decompressed")
    # error if the database doesn't exist
    if not Path(database_file).is_file():
        raise FileNotFoundError(database_file, "Database not found")
    
    qc_sql = Path("scripts/generate_reports.sql")
    if not qc_sql.is_file():
        raise FileNotFoundError(qc_sql, "generate_reports.sql QC SQL script not found")
    sql = qc_sql.read_text()

    con = duckdb.connect('output/monarch-kg.duckdb')
    con.execute(sql)

def export_tsv():
    export()


def do_prepare_release(dir: str = OUTPUT_DIR):

    compressed_artifacts = [
        'output/monarch-kg.duckdb',
        'output/monarch-kg-denormalized-edges.tsv',
        'output/monarch-kg-denormalized-nodes.tsv',
    ]

    for artifact in compressed_artifacts:
        if Path(artifact).exists() and not Path(f"{artifact}.gz").exists():
            sh.pigz(artifact, force=True)

    jsonl_tar = tarfile.open("output/monarch-kg.jsonl.tar.gz", "w:gz")
    jsonl_tar.add("output/monarch-kg_nodes.jsonl", arcname="monarch-kg_nodes.jsonl")
    jsonl_tar.add("output/monarch-kg_edges.jsonl", arcname="monarch-kg_edges.jsonl")
    jsonl_tar.close()

    os.remove("output/monarch-kg_nodes.jsonl")
    os.remove("output/monarch-kg_edges.jsonl")


def do_release(dir: str = OUTPUT_DIR, kghub: bool = False):

    # ensure that files that should be compressed are

    with open(f"{dir}/metadata.yaml", "r") as f:
        versions = yaml.load(f, Loader=yaml.FullLoader)

    release_ver = versions["kg-version"]

    logger = get_logger()
    logger.info(f"Creating dated release: {release_ver}...")

    try:
        logger.debug(f"Uploading release to Google bucket...")
        sh.touch(f"{dir}/{release_ver}")

        # copy data to monarch-archive bucket
        sh.gsutil(*f"-q -m cp -r data/* gs://monarch-archive/monarch-ingest-data-cache/{release_ver}".split(" "))

        # copy to monarch-archive bucket
        sh.gsutil(*f"-q -m cp -r {dir}/* gs://monarch-archive/monarch-kg-dev/{release_ver}".split(" "))

        # copy to data-public bucket
        sh.gsutil(
            *"-q -m cp -r".split(" "),
            f"gs://monarch-archive/monarch-kg-dev/{release_ver}",
            f"gs://data-public-monarchinitiative/monarch-kg-dev/{release_ver}",
        )

        # update "latest"
        sh.gsutil(*"-q -m rm -rf gs://data-public-monarchinitiative/monarch-kg-dev/latest".split(" "))
        sh.gsutil(
            *"-q -m cp -r".split(" "),
            f"gs://data-public-monarchinitiative/monarch-kg-dev/{release_ver}",
            "gs://data-public-monarchinitiative/monarch-kg-dev/latest",
        )

        # copy data to monarch-archive bucket only, so that we don't keep so many copies of the same huge files
        sh.gsutil(*f"-q -m cp -r data/* gs://monarch-archive/monarch-kg-dev/{release_ver}/data".split(" "))

        # index and upload to kghub s3 bucket
        if kghub:
            kghub_release_ver = str(release_ver).replace("-", "")
            logger.info(f"Uploading to kghub: {kghub_release_ver}...")

            # index files locally and upload to s3
            sh.multi_indexer(
                *f"-v --directory {dir} --prefix https://kghub.io/kg-monarch/{kghub_release_ver} -x -u".split(" ")
            )
            kg_hub_files = f"{dir}/monarch-kg.tar.gz {dir}/rdf/ {dir}/merged_graph_stats.yaml"
            sh.gsutil(
                *f"-q -m cp -r -a public-read {kg_hub_files} s3://kg-hub-public-data/kg-monarch/{kghub_release_ver}".split(
                    " "
                )
            )
            sh.gsutil(
                *f"-q -m cp -r -a public-read {kg_hub_files} s3://kg-hub-public-data/kg-monarch/current".split(" ")
            )
            # index files on s3 after upload
            sh.multi_indexer(
                *f"-v --prefix https://kghub.io/kg-monarch/ -b kg-hub-public-data -r kg-monarch -x".split(" ")
            )
            sh.gsutil(*f"-q -m cp -a public-read ./index.html s3://kg-hub-public-data/kg-monarch".split(" "))

        logger.debug("Cleaning up files...")
        sh.rm(f"output/{release_ver}")

        logger.info(f"Successfuly uploaded release!")
    except BaseException as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        logger.error(f"Oh no! Something went wrong:\n{fname}:{exc_tb.tb_lineno} - {exc_type} - {exc_obj}")
        logger.error(f"Traceback: \n{e.stderr}")

def clone_release(release_ver: str,
                  light: bool = False, 
                  dev: bool = False):

    logger = get_logger()
    logger.info(f"Cloning release: {release_ver}...")


    os.makedirs("output", exist_ok=True)
    os.makedirs("output/qc", exist_ok=True)

    if dev:
        kg_path = "monarch-kg-dev"
    else:
        kg_path = "monarch-kg"

    if light:        
        sh.gsutil(*f"-q -m cp -r gs://data-public-monarchinitiative/{kg_path}/{release_ver}/monarch-kg.tar.gz output/".split(" "))
        sh.gsutil(*f"-q -m cp -r gs://data-public-monarchinitiative/{kg_path}/{release_ver}/qc/monarch-kg-dangling-edges.tsv.gz output/".split(" "))

    else:
        sh.gsutil(*f"-q -m cp -r gs://data-public-monarchinitiative/{kg_path}/{release_ver}/* output/".split(" "))
