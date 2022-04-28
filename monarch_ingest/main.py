import subprocess, datetime
from kghub_downloader.download_utils import download_from_yaml
from monarch_ingest.cli_utils import *

import typer

typer_app = typer.Typer()

import logging
logging.basicConfig()
LOG = logging.getLogger(__name__)

OUTPUT_DIR = "output"

@typer_app.command()
def download(
    tags: Optional[List[str]] = typer.Option(None,  help="Which tags to download data for"),
    all: bool = typer.Option(False, help="Download all ingest datasets")
    ):
    if tags:
        download_from_yaml(
            yaml_file='monarch_ingest/download.yaml',
            output_dir='.',
            tags=tags,
        )
    elif all:
        download_from_yaml(
            yaml_file='monarch_ingest/download.yaml',
            output_dir='.'
        )

@typer_app.command()
def transform(
    output_dir: str = typer.Option(OUTPUT_DIR, help="Directory to output data"),
    # data_dir: str = typer.Option('data', help='Path to data to ingest),
    tag: str = typer.Option(
        None, help="Which ingest to run (see ingests.yaml for a list)"
    ),
    ontology: bool = typer.Option(False, help="Option: pass to run the ontology ingest"),
    all: bool = typer.Option(False, help="Ingest all sources"),
    row_limit: int = typer.Option(None, help="Number of rows to process"),
    do_merge: bool = typer.Option(False, "--merge", help="Merge output dir after ingest"),
    force: bool = typer.Option(
        None,
        help="Force ingest, even if output exists (on by default for single ingests)",
    ),
    quiet: bool = typer.Option(False, help="Suppress LOG output"),
    debug: bool = typer.Option(False, help="Print additional debug output to console"),
    log: bool = typer.Option(
        False, help="Write DEBUG level logs to ./logs/ for each ingest run"
    ),
):
    """
    Something descriptive
    """
    if ontology:
        LOG.info(f"Running ontology transform...")
        transform_ontology(output_dir, force)
    elif all:
        LOG.info(f"Running all ingests...")
        transform_all(
            f"{output_dir}/transform_output",
            row_limit=row_limit,
            force=force,
            quiet=quiet,
            debug=debug,
            log=log,
        )
    elif tag:
        if force is None:
            force = True
        transform_one(
            tag,
            f"{output_dir}/transform_output",
            row_limit=row_limit,
            force=force,
            quiet=quiet,
            debug=debug,
            log=log,
        )
    if do_merge:
        merge(f"{output_dir}/transform_output", output_dir)


@typer_app.command()
def merge(
    input_dir: str = typer.Option(
        f"{OUTPUT_DIR}/transform_output",
        help="Directory containing nodes and edges to be merged",
        ),
    output_dir: str = typer.Option(f"{OUTPUT_DIR}", help="Directory to output data"),
    ):
    """
    Something descriptive
    """
    merge_files(input_dir=input_dir, output_dir=output_dir)


@typer_app.command()
def release(
    update_latest: bool = typer.Option(False, help="Pass to update latest with this release")
):

    release_name = datetime.datetime.now()
    release_name = release_name.strftime("%Y-%m-%d")

    LOG.info(
        f"Creating release...\nToday's date: {release_name}"
    )

    try:
        LOG.debug(f"Uploading release to Google bucket...")
        subprocess.run(['touch', f"output/{release_name}"])
        subprocess.run(
            [
                'gsutil',
                '-m',
                'cp',
                '-r',
                'output/*',
                f"gs://monarch-ingest/{release_name}",
            ]
        )

        LOG.debug("Cleaning up files...")
        subprocess.run(['rm', f"output/{release_name}"])

        LOG.info(f"Successfuly uploaded release: see gs://monarch-ingest/{release_name}")
    except BaseException as e:
        LOG.error(f"Oh no! Something went wrong:\n{e}")

    if update_latest:
        LOG.debug(f"Replacing latest with this release")
        subprocess.run(['gsutil', '-m', 'rm', '-rf', 'gs://monarch-ingest/latest'])
        subprocess.run(['gsutil','-m','cp','-r',f"gs://monarch-ingest/{release_name}","gs://monarch-ingest/latest",])
        LOG.info(f"Updated 'latest' to current release.")

if __name__ == "__main__":
    typer_app()
