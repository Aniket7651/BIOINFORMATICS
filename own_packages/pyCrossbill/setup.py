import setuptools
  
with open("README.txt", "r") as file:
    description = file.read()
  
setuptools.setup(
    name="pyCrossbill",
    version="0.1",
    author="aniket_yadav",
    author_email="aniketyadav8687@gmail.com",
    packages=["pyCrossbill"],
    description="A bioinformatics tool for performing various task related to database, neucleotide, protein, and some compounds",
    long_description=description,
    long_description_content_type="text",
    url="https://github.com/gituser/test-tackage",
    license='Apache',
    python_requires='>=3.8',
    install_requires=[]
)