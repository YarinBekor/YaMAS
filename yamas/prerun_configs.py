from .utilities import run_cmd

command_based_on_os = {'Ubuntu': ['wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz'],
                       'CentOS': ['wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz'],
                       }


def set_environment(operating_system_type):
    command = command_based_on_os[operating_system_type]
    print("######   Downloading SRA Toolkit   ######")
    run_cmd(command)

    print("######   Extracting the contents of the tar file   ######")
    run_cmd(['tar -vxzf sratoolkit.tar.gz > sra_toolkit_download_log.txt'])
    with open('sra_toolkit_download_log.txt', 'r') as file:
        last_line = file.readlines()[-1]
        version = last_line.split("/")[0]
        print(f"######   export PATH   ######")
        try:
            run_cmd([f'export PATH=$PATH:$PWD/{version}/bin'])
        except Exception as e:
            print("######   Can't export the path. Failed. See below.   ######")

    print("######   Checking if we are ready to go...   ######")
    run_cmd(['which fastq-dump > check_fastq-dump.txt'])
    with open('check_fastq-dump.txt', 'r') as check_file:
        content= check_file.read()
        if f'{version}/bin/fastq-dump' in content:
            print("######   YAMAS is ready to run!   ######")
        else:
            print("######   YAMAS is NOT ready! Try again or you can set the environment by yourself.\n  Follow the instructions (from step 2) that can be found in git: https://github.com/YarinBekor/YaMAS.   ######\n    ")
