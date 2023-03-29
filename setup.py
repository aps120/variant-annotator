from setuptools import setup

# read requirements.
requirements = []
with open("requirements.txt", 'rU') as reader:
    for line in reader:
        requirements.append(line.strip())

setup(name='vcfAnnotator',
        python_requires='>=3',
        version='20230327',
        description='Annotate vcf files',
        license='Apache 2.0',
        packages=['vcfAnnotator'],
        install_requires=requirements,
        entry_points = {
            'console_scripts' : ['vcf-annotator=vcfAnnotator.vcf_annotator:main'],
            },
        zip_safe=False)
