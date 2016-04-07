import setuptools
from os import path



whereAmI = os.path.dirname(os.path.realpath(__file__))




setup(
name='rose2'
version='1.0.0'
description='ROSE2 python package'
long_description='PROGRAM TO STITCH TOGETHER REGIONS TO FORM ENHANCERS, MAP READ DENSITY TO STITCHED REGIONS,AND RANK ENHANCERS BY READ DENSITY TO DISCOVER SUPER-ENHANCERS'
url=''
author=''
license=''

classifiers=[

	'Development Status :: 3 - Alpha',

	]
keywords='stitched regions enhancer reads',

install_requires=[''],

extras_require={},


entry_points={
    'console_scripts': [
        'rose2=rose2:rose2',
        ]
    }

)