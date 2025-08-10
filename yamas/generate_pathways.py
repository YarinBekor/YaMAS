import os
from pathlib import Path
from typing import Union, List, Optional
from .utilities import run_cmd


def run_humann_pipeline(
    input_file: Union[str, Path],
    output_dir: Union[str, Path],
    meta_profile: Union[str, Path],
    threads: int = 8,
    chocophlan_db: Optional[Union[str, Path]] = None,
    uniref_db: Optional[Union[str, Path]] = None,
    utility_db: Optional[Union[str, Path]] = None,
    resume: bool = False,
    input_format: str = "fastq"
) -> List[Path]:
    """
    Run the HUMAnN pipeline from the middle, using existing MetaPhlAn and Bowtie outputs.
    Supports FASTQ inputs or compressed SAM (.sam.bz2) with bypass of Bowtie2 steps.
    """
    input_path = Path(input_file)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[HUMAnN] Starting pipeline for {input_path.name}")
    print(f"[HUMAnN] Input format: {input_format}")
    
    # Build base HUMAnN command
    cmd_parts = [
        "humann",
        f"--output {out_dir}",
        f"--threads {threads}",
        f"--taxonomic-profile {meta_profile}",
        f"--input-format {input_format}"
    ]

    # Adjust input and bypass flags
    if input_format == "fastq":
        print(f"[HUMAnN] Using FASTQ input: {input_path.name}")
        cmd_parts.append(f"--input {input_path}")
    
    elif input_format == "sam":
        # only treat true .sam.bz2, else fallback to fastq
        if input_path.name.endswith(".sam.bz2"):
            sam_path = input_path.with_suffix("").with_suffix(".sam")
            if not sam_path.exists():
                print(f"[HUMAnN] Decompressing {input_path.name} to {sam_path.name}")
                run_cmd([f"bzip2 -dc {input_path} > {sam_path}"])
                print(f"[HUMAnN] Decompression complete: {sam_path.name}")
            cmd_parts.append(f"--input {sam_path}")
            cmd_parts.extend(["--bypass-nucleotide-index", "--bypass-prescreen", "--bypass-nucleotide-search"])
            print(f"[HUMAnN] Prepared SAM bypass on: {sam_path.name}")
        else:
            # fallback to fastq branch
            input_format = "fastq"
            cmd_parts[cmd_parts.index(f"--input-format {input_format}")] = f"--input-format fastq"
            cmd_parts.append(f"--input {input_path}")
    else:
        print(f"[HUMAnN] Using other input: {input_path.name}")
        cmd_parts.append(f"--input {input_path}")

    # Optional databases
    if chocophlan_db:
        cmd_parts.append(f"--nucleotide-database {chocophlan_db}")
    if uniref_db:
        cmd_parts.append(f"--protein-database {uniref_db}")
    if utility_db:
        cmd_parts.append(f"--utility-map {utility_db}")
    if resume:
        cmd_parts.append("--resume")

    # Prepare log redirection
    log_file = out_dir / "humann.log"
    full_cmd = " ".join(cmd_parts) + f" > {log_file} 2>&1"
    print(f"[HUMAnN] Running command: {full_cmd}")

    # Execute via run_cmd
    run_cmd([full_cmd])
    print(f"[HUMAnN] Finished command for {input_path.name}")

    # Return and report pathways output files
    outputs = sorted(out_dir.glob("*_pathabundance.tsv"))
    print(f"[HUMAnN] Found {len(outputs)} pathway files: {[p.name for p in outputs]}")
    return outputs
