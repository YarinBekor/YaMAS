from utilities import run_cmd

command_based_on_os = {'Ubuntu': ['wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz'],
                       'CentOS': ['wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz'],
                       'MacOS': ['curl --output sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz']
                       }


def set_enviorment(operating_system_type):
    command = command_based_on_os[operating_system_type]
    run_cmd(command)
    run_cmd(['tar -vxzf sratoolkit.tar.gz > sra_toolkit_download_log.txt'])
    with open('sra_toolkit_download_log.txt', 'r') as file:
        last_line = file.readlines()[-1]
        version = last_line.split("/")[0]
        run_cmd([f'export PATH=$PATH:$PWD/{version}/bin'])
    
    print("######   YAMAS is ready to run!   ######")
