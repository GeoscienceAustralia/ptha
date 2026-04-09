# ptha
Codes for Probabilistic Tsunami Hazard Assessment (PTHA). 

The folder [ptha_access](ptha_access) contains scripts/tutorials to access the
[2018 Australian PTHA database](http://dx.doi.org/10.11636/Record.2018.041).

The folder [misc](misc) contains data and code from GA's related tsunami
projects, with links to the associated papers and technical reports.

The folder [propagation](propagation) contains the nested grid shallow water
equations solver [SWALS](./propagation/SWALS/), used to model tsunamis in
many of our studies.

The folder [R](R) includes the R package [rptha](R/rptha) along with
[installation instructions](R/README.md). It also contains various 
[tutorials and template scripts](R/examples) that use
rptha, (including 
[project-specific scripts used for the 2018 Australian PTHA](R/examples/austptha_template/) ).

*Please use the most up-to-date source code, not the "releases".  In this
repository, "releases" are only used to snapshot the code at some point in time
(e.g. a reference for older studies). The master branch is recommended for general
use, and should pass the full test suite.*


## License

This code is licensed under a BSD 3-clause license. See the [license deed](LICENSE).

## Contacts

**Gareth Davies**  
*Lead Developer*  
Geoscience Australia
<gareth.davies@ga.gov.au>
