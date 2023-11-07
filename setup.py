import setuptools
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setuptools.setup(
    name="YMS",
    long_description=long_description,
    long_description_content_type='text/markdown',
    version="1.1.4",
    author="Yarin Bekor",
    author_email="yarin.bekor@gmail.com",
    description="YOLO Microbiome Analysis System",
    license='MIT',
    entry_points={
        'console_scripts': [
            'yamas = yamas:main',
        ]
    },
    packages=setuptools.find_packages(),
    package_data={
        'yamas': ['config.json']
    },
    install_requires=[
        'tqdm',
        'MetaPhlAn'
    ]
)