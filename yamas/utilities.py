import os
import subprocess
from dataclasses import dataclass


@dataclass(frozen=True)
class ReadsData:
    dir_path: str
    fwd: bool = True
    rev: bool = False


def run_cmd(command: list):
    os.system(" ".join(command))


def qiime2_version():
    o, e = run_cmd(["conda", "env", "list"])
    qiime_version = ""
    for env_1 in o.split("\n"):
        for env_2 in env_1.split("/"):
            if "qiime2" in env_2 and " " not in env_2:
                qiime_version = env_2
    return qiime_version


def download_classifier_url():
    return f"https://data.qiime2.org/{qiime2_version().split('-')[1]}/common/gg-13-8-99-nb-classifier.qza"


def check_conda_qiime2():
    conda_prefix = os.environ.get("CONDA_PREFIX", None)
    if conda_prefix is None:
        raise Exception("Tha package can work only on CONDA environment QIIME2.")
    if os.path.split(conda_prefix)[-1].split("-")[0] != "qiime2":
        raise EnvironmentError("The package must work on qiime2 environment. "
                               "\nActivate the requested environment, "
                               "or create it with 'conda install -c qiime2 qiime2'")