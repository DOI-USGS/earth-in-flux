# USGS Land Change Chart Gallery

> _A newer version of the software may be available. See https://code.usgs.gov/wma/vizlab/earth-in-flux/-/releases to view all releases._

The land change chart gallery is a shared initiative between the USGS Water Mission Area and Ecosystems Mission Area to communicate key findings of USGS land change science in innovative ways, and to encourage creativity, exploration, and community in data visualization across USGS mission areas.

**The chart gallery can be viewed at [https://labs.waterdata.usgs.gov/visualizations/earth-in-flux](https://labs.waterdata.usgs.gov/visualizations/earth-in-flux).**

## Citation

Corson-Dosch, Hayley, Archer, Althea, Kwang, Jeffrey, Jaenicke, Margaret, and Nell, Cee. 2024. USGS Land Change Chart Gallery. U.S. Geological Survey software release. Reston, VA. [https://doi.org/10.5066/P1J2PMFE](https://doi.org/10.5066/P1J2PMFE)

## Additional contributors
[Laura Gemery](https://www.usgs.gov/staff-profiles/laura-gemery) and [Natalie Kerhwald](https://www.usgs.gov/staff-profiles/natalie-m-kehrwald) consulted on the development of this website as subject matter experts.

## Additional information
* We welcome contributions from the community. See the [guidelines for contributing](https://github.com/DOI-USGS/earth-in-flux/) to this repository on GitHub.
* [Disclaimer](https://code.usgs.gov/wma/vizlab/earth-in-flux/-/blob/main/DISCLAIMER.md)
* [License](https://code.usgs.gov/wma/vizlab/earth-in-flux/-/blob/main/LICENSE.md)

## To build the website locally
Clone the repo. In the directory, run `npm install` to install the required modules. This repository requires `npm v20` to run. If you are using a later version of `npm`, you may [try using `nvm` to manage multiple versions of npm](https://betterprogramming.pub/how-to-change-node-js-version-between-projects-using-nvm-3ad2416bda7e).

Once the dependencies have been installed, run `npm run dev` to run the site locally from your browser.

## Notes for adding content to the site
The master file controlling viz content is `'src/assets/content/ChartGrid.js'`. This controls the cards seen in the landing view, as well as the content used to populate the visualization subpages. Note that the page routing setup requires that very specific naming conventions be followed. In `'src/assets/content/ChartGrid.js'`, there is an object for each visualization. The critical parameter to note is `vizKey`, which is used to dynamically load the component, text, references, and authors content for each visualization subpage. The component associated with the visualization must be named `src/components/${vizKey}Viz.vue` (e.g., `'src/components/GlacierScanViz.vue'`), and the nested objects in `'src/assets/text/text.js'` and `'src/assets/text/authors.js'` (see notes, below), must be named to match `vizKey`.

### Page text setup
  * All page text should be added directly to `'src/assets/text/text.js'`. The master text object from that file is automatically imported in `'src/views/VisualizationView.vue'` and used to set the page title. The text object is also imported in `'src/components/SubPage.vue'` and nested objects containing text for each visualization page are dynamically passed to visualization subpages. The nested objects must be named to match the `vizKey` specified for each viz in `ChartGrid.js`

### Authors setup
  * All authors for your viz should be added directly to `'src/assets/text/authors.js'`. There is a object for each visualization. The content will be dynamically passed to visualization subpages. The objects must be named to match the `vizKey` specified for each viz in `ChartGrid.js`.

### References setup
  * All references for a project should be added directly to `'src/assets/text/references.js'`. There is a object for each project. The content will be dynamically passed to project subpages.

### Data Pipelines
Any preparatory pipeline code should be placed in subdirectories in the main repo. Please include a README.md in that subdirectory with instructions on how to run your code.
