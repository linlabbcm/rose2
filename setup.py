import os

from setuptools import setup, find_packages



setup(
	name='rose2',
	version='1.0.0',
	description='ROSE2 python package',
	long_description='PROGRAM TO STITCH TOGETHER REGIONS TO FORM ENHANCERS, MAP READ DENSITY TO STITCHED REGIONS,AND RANK ENHANCERS BY READ DENSITY TO DISCOVER SUPER-ENHANCERS',
	url='',
	author='',
	license='',

	classifiers=[

		],
	keywords='stitched regions enhancer reads',

	packages=find_packages(),
	include_package_data=True,
	package_data={},

	install_requires=[],

	extras_require={},


	scripts=['scripts/ROSE2_callSuper.R', 'scripts/ROSE2_stitchOpt.R'],

	entry_points={
	    'console_scripts': [
	        'rose2=rose2:main',
	        ]
	    },

)