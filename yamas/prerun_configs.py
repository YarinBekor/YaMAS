from .utilities import run_cmd

command_based_on_os = {'Ubuntu': ['wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz'],
                       'CentOS': ['wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz'],
                       }


def set_environment(operating_system_type):

    command = command_based_on_os[operating_system_type]
    try:
        run_cmd(command)
    except Exception as e:
        print(f"######  Problem with downloading SRA Toolkit. \nCheck if this is your correct operating system. You typed: {operating_system_type}.Try again or you can set the environment by yourself.\nFollow the instructions (from step 2) that can be found in git: https://github.com/YarinBekor/YaMAS.   ######\n {e}")


    try:
        run_cmd(['tar -vxzf sratoolkit.tar.gz > sra_toolkit_download_log.txt'])
        with open('sra_toolkit_download_log.txt', 'r') as file:
            last_line = file.readlines()[-1]
            version = last_line.split("/")[0]
            run_cmd([f'export PATH=$PATH:$PWD/{version}/bin'])
    except Exception as e:
        print(f"######   Problem with extracting the content of the tar file. \nTry again or you can set the environment by yourself.\n  Follow the instructions (from step 2) that can be found in git: https://github.com/YarinBekor/YaMAS.   ######\n{e} ")


        run_cmd(['which fastq-dump > check_fastq-dump.txt'])
        with open('check_fastq-dump.txt', 'r') as file:
            content= file.read()
            if f'{version}/bin/fastq-dump' in content:
                print("######   YAMAS is ready to run!   ######")
            else:
                print("######   YAMAS is NOT ready! Try again or you can set the environment by yourself.\n  Follow the instructions (from step 2) that can be found in git: https://github.com/YarinBekor/YaMAS.   ######\n{e}    ")
