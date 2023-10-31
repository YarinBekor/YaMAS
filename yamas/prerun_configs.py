from .utilities import run_cmd

command_based_on_os = {'Ubuntu': ['wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz'],
                       'CentOS': ['wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz'],
                       }


def set_enviorment(operating_system_type):
    command = command_based_on_os[operating_system_type]
    run_cmd(command)
    with open('sra_toolkit_download_log.txt', 'w') as file:
        pass
    run_cmd(['tar -vxzf sratoolkit.tar.gz > sra_toolkit_download_log.txt 2>&1'])
    with open('sra_toolkit_download_log.txt', 'r') as file:
        last_line = None
        for line in file:
            last_line = line
        version = last_line.split("/")[0]
        run_cmd([f'export PATH=$PATH:$PWD/{version}/bin'])
    
    print("######   YAMAS is ready to run!   ######")
