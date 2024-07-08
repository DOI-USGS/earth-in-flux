# Climate Chart Gallery Contribution

Please fill out this information prior to creating a merge request, *unless you've filled out the template R or Python scripts*. After your merge request is approved, please make sure the image linked here is your final image in the `out/` folder (the `out/` png file itself is not committed, but the README with the linked image is). We will also ask you to share the final image in a sharepoint folder provided via gitlab.

Your name: Althea Archer

Your project: Beaufort Sea

Image:

Share the link in this document with the format `![](out/BeaufortSeaTimeline.png)`:

![](out/BeaufortSeaTimeline.png)

Caption with key message communicated by this viz :

1. 2000 years ago, the system was dominated by Paracyprideis, Cassidulina and Elphidium species.
2. The first 1000 years demonstrated relative stability
3. Then, from 1000-1300, agglutinated species Spiroplectammina became more abundant as a result of the climate instability during the Medieval Climate Anomaly
4. During the Little Ice Age, Paracyprideis declined to its lowest levels in the record
5. In the modern record, Spiroplectammina and Kotoracythere are at all-time highs, and the previously dominant Elphidium and Cassidulina are sparser

Your data sources: 

Sciencebase:

> Gemery, L., 2021, Data Release to Multi-proxy record of ocean-climate variability during the last 2 millennia on the Mackenzie Shelf, Beaufort Sea (2013): U.S. Geological Survey data release, https://doi.org/10.5066/P9SRRW6T.

Original manuscript: 

> Gemery, L., Cronin, T.M., Cooper, L.W., Roberts, L.R., Keigwin, L.D., Addison, J.A., Leng, M.J., Lin, P., Magen, C., Marot, M.E., and Schwartz, V., 2023, Multi-proxy record of ocean-climate variability during the last two millennia on the Mackenzie Shelf, Beaufort Sea: Micropaleontology, v. 69, no. 3, p. 345–360, https://doi.org/10.47894/mpal.69.3.04.

Foram data: 

> Seidenstein, J.L., Cronin, T.M., Gemery, L. et al. Late Holocene paleoceanography in the Chukchi and Beaufort Seas, Arctic Ocean, based on benthic foraminifera and ostracodes. Arktos 4, 1–17 (2018). https://doi.org/10.1007/s41063-018-0058-7


Overall method to create this viz:

1. Create an RStudio project in the `beaufortSea` folder
2. Download data sources and store them in the `beaufortSea/in` folder
3. Run `targets::tar_make()` to create the visualization
