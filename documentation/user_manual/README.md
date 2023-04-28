# BerkeleyGW user manual

This directory contains the BerkeleyGW documentation for the end user.

Please refer to the [BerkeleyGW main page](https://berkeleygw.org/) for the
full documentation if you would like to skip making the manual yourself.


## External documentation

Some external links for resources, including the online version of the
BerkeleyGW manual:

- [BerkeleyGW manual](http://manual.berkeleygw.org/)
- [BerkeleyGW tutorial](https://berkeleygw.org/tutorial/)
- [BerkeleyGW implementation paper](https://arxiv.org/abs/1111.4429)

Previous workshops:

- [2019](https://sites.google.com/site/berkeleygw2019/)
- [2018](https://sites.google.com/site/berkeleygw2018/)
- [2017](https://sites.google.com/site/berkeleygw2017/home)
- [2015](https://sites.google.com/site/berkeleygw2015/home)
- [2014](https://sites.google.com/site/berkeleygw2014/home)


## Building the manual locally

Before you build the manual locally in your workstation, you need to have
`mkdocs` installed, which can achieved with:

    pip install --user mkdocs mkdocs-material python-markdown-math mako

You will need `mkdocs>=1.0`; we tested the following version of the packages:

    mkdocs                        1.0.4
    mkdocs-material               4.0.2
    python-markdown-math          0.3

To build the manual, type:

    python assemble_manual.py


To see the manual, go into the `website/` directory and type:

    mkdocs serve


You can then open a browser and go to the IP address indicated by `mkdocs`,
which is probably:

    http://127.0.0.1:8000


#BEGIN_INTERNAL_ONLY
## Publishing manual (developers only)

To publish the manual, go to the website directory and issue:

    mkdocs build
    aws s3 cp site/ s3://manual.berkeleygw.org/2.2/ --recursive

You will need to have preconfigured `aws` in your machine, of course, with the
correct access keys. You can also completely remove any previous manual with:

    aws s3 rm s3://manual.berkeleygw.org/2.2/ --recursive

Change the `2.2` above as needed.
#END_INTERNAL_ONLY
