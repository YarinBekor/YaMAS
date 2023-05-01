import setuptools

setuptools.setup(
    name="YMS",
    version="0.71",
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
        'tqdm'
    ]
)