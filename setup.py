import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name='TAP',
	version='v0.0.1',
	author="Kenny Pavan",
	author_email="pavan@ohsu.edu",
	description="A Python tool for feature selection of multiplex single-cell tropism experiments. ",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/ArpiarSaundersLab/tap",
	packages=setuptools.find_packages(where='src'),  
	package_dir={'': 'src'}, 
    include_package_data=True,
    package_data={
        "tap": ["templates/*.html","models/*"]
    },
	python_requires='>=3.12',
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	install_requires=[
		'scanpy>=1.10.3',
		'annsql>=1.0.0',
	],
)
