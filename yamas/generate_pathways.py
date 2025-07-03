import subprocess
from typing import Union, List, Optional
from pathlib import Path

def run_humann_pipeline(fastq_input: Union[str, Path],
                        output_dir: Union[str, Path],
                        meta_profile: Union[str, Path],
                        threads: int = 8,
                        chocophlan_db: Optional[Union[str, Path]] = None,
                        uniref_db: Optional[Union[str, Path]] = None,
                        utility_db: Optional[Union[str, Path]] = None,
                        resume: bool = False
                        ) -> List[Path] :
    '''
    Run HUMAnN pipeline on the given FASTQ file.
    
    Parameters:
    - fastq_path: Path to the input FASTQ file
    - output_dir: Directory where HUMAnN output will be saved
    - threads: Number of threads to use (default: 4)
    - chocophlan_db: Path to the ChocoPhlAn database (optional)
    - uniref_db: Path to the UniRef database (optional)
    - utility_db: Path to the utility database (optional)
    - resume: Whether to resume from a previous run (default: False)
      Returns
    -------
    List[Path]
        List of result files produced by HUMAnN.
        
    '''
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "humann",
        "--input", str(fastq_input),
        "--output", str(output_dir),
        "--threads", str(threads)
    ]
    if meta_profile:
        cmd += ["--taxonomic-profile", str(meta_profile)]
    
    # Add the database paths if provided
    if chocophlan_db:
        cmd += ["--nucleotide-database", str(chocophlan_db)]
    if uniref_db:
        cmd += ["--protein-database", str(uniref_db)]
    if utility_db:
        cmd += ["--utility-map", str(utility_db)]
    if resume:
        cmd += ["--resume"]
    
    log_file = output_dir / "humann.log"
    with log_file.open("w") as log_fp:
        result = subprocess.run(cmd, stdout=log_fp,
                                stderr=subprocess.STDOUT, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"HUMAnN failed - see log at {log_file}")

    print("HUMAnN pipeline completed successfully.")
    tables = sorted(output_dir.glob("*_pathabundance.tsv"))
    
    return tables