import setuptools

setuptools.setup(
    name="YMS",
    version="0.1",
    author="Yarin Bekor",
    author_email="yarin.bekor@domain.com",
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
    }
)
