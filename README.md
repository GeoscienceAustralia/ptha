# ptha
Codes for probabilistic tsunami hazard assessment. 

The folder [ptha_access](ptha_access) contains scripts/tutorials to access the
2018 Australian PTHA results - or typically subsets thereof, since the full
analysis is several TB in size.

The folder [misc](misc) contains data and code for related work, as detailed therein.

The folder [propagation](propagation) contains a shallow water equations solver.

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
