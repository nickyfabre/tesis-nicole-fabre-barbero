<div class="cover">

<p class="logo"><img src="logo_fcen.png" alt="FCEN-UBA" width="160"></p>

<p class="inst">Universidad de Buenos Aires<br>Facultad de Ciencias Exactas y Naturales</p>

<p class="kind">Tesis para optar por el título de Licenciada en Ciencias Biológicas</p>

<p class="title">Asignación filogenética proporcional de péptidos compartidos para la estimación de biomasa microbiana en metaproteómica</p>

<p class="people"><strong>Autora:</strong> Nicole Fabre Barbero<br><strong>Director:</strong> Dr. Lukas Käll<br><strong>Codirector:</strong> Dr. Germán Rosano</p>

<p class="place"><strong>Lugar de trabajo:</strong> Science for Life Laboratory (SciLifeLab), KTH Royal Institute of Technology, Estocolmo, Suecia (trabajo remoto)</p>

<p class="date">Octubre de 2026</p>


</div>

# Resumen {.unnumbered}

La metaproteómica permite estimar la contribución de distintas especies a la biomasa total de una comunidad microbiana, pero dos prácticas habituales limitan estas estimaciones. Los péptidos compartidos entre organismos cercanos se descartan o se asignan sólo al ancestro común más cercano (LCA, por sus siglas en inglés) de esos organismos, y las identificaciones que no superan un umbral fijo de confianza se eliminan, mientras que las restantes se tratan como correctas. Esta tesis se pregunta si los péptidos compartidos pueden usarse para estimar la abundancia de cada especie, y si la probabilidad de que cada identificación sea correcta puede reemplazar al umbral de confianza. Se compararon cinco estrategias que asignan las intensidades de los péptidos sobre un árbol de especies, cada una de las cuales agrega un elemento a la anterior: descartar los péptidos compartidos, conservarlos en el LCA, repartirlos entre las cepas candidatas una vez o de forma iterativa con un algoritmo de esperanza-maximización (EM) y, por último, ponderar cada identificación por su probabilidad de ser correcta, 1 − PEP, donde PEP es la probabilidad de error posterior. Las estrategias se evaluaron en comunidades artificiales de unas 30 cepas de composición conocida, procesadas sin filtro de tasa de falsos descubrimientos, y en comunidades simuladas en las que se variaron la proporción de péptidos compartidos y la tasa de identificaciones falsas. Conservar los péptidos compartidos en el LCA dio totales correctos para los clados, mientras que descartarlos subestimó fuertemente los clados formados por cepas emparentadas. A nivel de especie, solo el reparto de la señal compartida corrigió las cepas con señal mayormente compartida, cuyo error pasó de decenas de veces a menos de tres veces. Usar todas las identificaciones ponderadas por 1 − PEP fue significativamente más exacto que aplicar un umbral estricto en dos de los tres diseños de la comunidad, principalmente porque recuperó las cepas poco abundantes, a costa de un leve aumento de la señal falsa en las cepas ausentes. En las simulaciones, cuando la confianza de las identificaciones correctas dependía de su intensidad, como en los datos reales, esta ponderación redujo el error de las especies raras en casi todas las comunidades, mientras que usar todas las identificaciones sin ponderarlas aumentó la señal falsa. Estos resultados muestran que los péptidos compartidos deben conservarse en el árbol y repartirse entre las especies, y que la confianza de cada identificación puede reemplazar a un umbral fijo al estimar abundancias.

**Palabras clave:** metaproteómica, péptidos compartidos, ancestro común más cercano, probabilidad de error posterior, estimación de biomasa, comunidad artificial.

# Abstract {.unnumbered}

**Proportional phylogenetic allocation of shared peptides for microbial biomass estimation in metaproteomics**

Metaproteomics allows for the estimation of the contribution of different species to the total biomass of a microbial community, but two common practices limit these estimates. Peptides shared by closely related organisms are either discarded or assigned only to the lowest common ancestor (LCA) of those organisms, and identifications below a fixed confidence threshold are removed, while the remaining ones are treated as correct. This thesis asks whether shared peptides can be used to estimate the abundance of each species, and whether the probability that each identification is correct can replace the confidence threshold. Five strategies that allocate peptide intensities on a species tree were compared, each adding one element to the previous one: discarding shared peptides, keeping them at the LCA, splitting them among the candidate strains once or iteratively with an expectation–maximization (EM) algorithm and, finally, weighting every identification by its probability of being correct, 1 − PEP, where PEP is the posterior error probability. The strategies were evaluated on mock communities of about 30 strains with known composition, processed without false discovery rate filtering, and on simulated communities in which the proportion of shared peptides and the rate of false identifications were varied. Keeping shared peptides at the LCA gave correct totals for the clades, whereas discarding them strongly underestimated the clades formed by closely related strains. At the species level, only splitting the shared signal corrected the strains whose signal is mostly shared, reducing their error from more than tenfold to less than threefold. Using every identification weighted by 1 − PEP was significantly more accurate than a strict threshold in two of the three community designs, mainly because it recovered low-abundance strains, at the cost of a small increase in false signal on absent strains. In the simulations, when the confidence of correct identifications depended on their intensity, as in the real data, this weighting reduced the error of rare species in almost every community, whereas using every identification without weights increased the false signal. These results show that shared peptides should be kept on the tree and split among species, and that the confidence of each identification can replace a fixed threshold when estimating abundances.

**Keywords:** metaproteomics, shared peptides, lowest common ancestor, posterior error probability, biomass estimation, mock community.
