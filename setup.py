import setuptools

setuptools.setup(
    name="YaMaS",
    version="0.1",
    author="Yarin Bekor",
    author_email="yarin.bekor@domain.com",
    description="YOLO lab Microbiome System",
    license='MIT',
    entry_points={
        'console_scripts': [
            'yamas = yamas:main',
        ]
    }
)