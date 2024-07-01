# climate-charts

The climate chart gallery is a shared initiative between the USGS Water Mission Area, Ecosystems Mission Area, and Climate Adaptation Science Centers to communicate key findings of USGS climate science in innovative ways, and to encourage creativity, exploration, and community in data visualization across USGS mission areas.

## To build the website locally
Clone the repo. In the directory, run `npm install` to install the required modules. This repository requires `npm v20` to run. If you are using a later version of `npm`, you may [try using `nvm` to manage multiple versions of npm](https://betterprogramming.pub/how-to-change-node-js-version-between-projects-using-nvm-3ad2416bda7e).

Once the dependencies have been installed, run `npm run dev` to run the site locally from your browser.

## Notes for adding content to the site/editing draft content
The master file controlling viz content is `'src/assets/content/ChartGrid.js'`. This controls the cards seen in the landing view, as well as the content used to populate the visualization subpages. Note that the page routing setup requires that very specific naming conventions be followed. In `'src/assets/content/ChartGrid.js'`, there is an object for each visualization. The critical parameter to note is `vizKey`, which is used to dynamically load the component, text, references, and authors content for each visualization subpage. The component associated with the visualization must be named `src/components/${vizKey}Viz.vue` (e.g., `'src/components/GlacierScanViz.vue'`), and the nested objects in `'src/assets/text/text.js'`, `'src/assets/text/references.js'`, and `'src/assets/text/authors.js'` (see notes, below), must be named to match `vizKey`.

Page text setup
  * All page text should be added directly to `'src/assets/text/text.js'`. The master text object from that file is automatically imported in `'src/views/VisualizationView.vue'` and used to set the page title. The text object is also imported in `'src/components/SubPage.vue'` and nested objects containing text for each visualization page are dynamically passed to visualization subpages. The nested objects must be named to match the `vizKey` specified for each viz in `ChartGrid.js`

References setup
  * All references for your viz should be added directly to `'src/assets/text/references.js'`. There is a object for each visualization. The content will be dynamically passed to visualization subpages. The objects must be named to match the `vizKey` specified for each viz in `ChartGrid.js`. See the `TEMPLATE` object for how to format different types of references.

Authors setup
  * All references for your viz should be added directly to `'src/assets/text/authors.js'`. There is a object for each visualization. The content will be dynamically passed to visualization subpages. The objects must be named to match the `vizKey` specified for each viz in `ChartGrid.js`. See the `TEMPLATE` object for how to specify different authors. Note that we may need to refine how we list authors/describe their contributions.
 
## General notes for development when using this template
1. This website template uses Vue 3 and the `<script setup>` composition API syntax to build components, which requires less boilerplate. See the [`<script setup>` guide](https://vuejs.org/api/sfc-script-setup.html). Any top-level defined variables or imported components are directly available for use in the `<template>`. Components now no longer need to be explicitly named, and can be imported directly by name using the filename, e.g. `import HeaderUSWDSBanner from "@/components/HeaderUSWDSBanner.vue"`.
2. Please do not delete or make any modifications to the following components/files:
    * The `node_modules` directory
    * `'public/favicon.ico'`
    * In `'src/assets'`:
      * The USGS Viz ID Stylesheets in `'src/assets/css'`:
        * `'common.css'`
        * `'custom.css'` 
      * The contents of the `'src/assets/usgsHeaderAndFooter'` directory
    * In `'src/components/'`
        * `'AuthorshipSection.vue'` - this should only be edited by the webdev team
        * `'ChartCard.vue'` - this should only be edited by the webdev team
        * `'ChartGrid.vue'` - this should only be edited by the webdev team
        * `'FooterUSGS.vue'`
        * `'HeaderUSGS.vue'`
        * `'HeaderUSWDSBanner.vue'`
        * `'NavBar.vue'` - this should only be edited by the webdev team
        * `'PreFooterCodeLinks.vue'`
        * `'ReferencesSection.vue'` - this should only be edited by the webdev team
        * `'SubPage.vue'` - this should only be edited by the webdev team
        * `'VizSection.vue'` _(see caveat in 4., below)_
        * `'WindowSize.vue'`
        * `'WorkInProgressWarning.vue'`
    * `'src/router/index.js'`
    * `'src/stores/WindowSizeStore.js'`
    * `'src/views/Error404Page.vue'`
    * `'src/App.vue'`
    * `'.env.beta_tier'`
    * `'.env.prod_tier'`
    * `'.env.test_tier'`
    * `'.eslintrc.cjs'`
    * `'.prettierrc.json'`
    * `'build.sh'`
    * `'CODE_OF_CONDUCT.md'`
    * `'DISCLAIMER.md'`
    * `'Dockerfile'`
    * `'LICENSE.md'`
    * `'package-lock.json'` _Note: will be modified indirectly by edits to `'package.json'`, and should be committed and pushed, but should not be edited directly_
3. The following files should have only limited edits, as specified in [Steps when using as template for new project](#steps-when-using-as-template-for-new-project), above:
    * `'jenkins/Jenkinsfile.build'`
    * In `'src/assets/text'`:
      * `'authors.js'`
      * `'references.js'`
    * `'.env'`
    * `'code.json'`
    * `'contributing.md'`
    * `'index.html'` (will also need to edit to modify analytics setup, per the [Vizlab website release checklist](https://doimspp.sharepoint.com/:w:/r/sites/IIDDStaff/_layouts/15/Doc2.aspx?action=edit&sourcedoc=%7B3c0899c4-cc87-4c82-a7e2-3f8e78439083%7D&wdOrigin=TEAMS-MAGLEV.teamsSdk_ns.rwc&wdExp=TEAMS-TREATMENT&wdhostclicktime=1714053079214&web=1))
    * `'package.json'` (may need additional edits if more packages are needed, see 3., below) 
4. Depending on what visualization tools you are using, you may need to add additional packages, which will require edits to some or all of the following files (reach out to Hayley with questions):
    * `'package.json'`
      * Any edits to `'package.json'` will require re-running `npm install`, which will update `'package-lock.json'`. **Be sure to commit and push any changes to `'package-lock.json'`**.
    * `'src/main.js'` 
    * `'vite.config.mjs'`. 
5. The template component `'src/components/VizSection.vue'` is designed to be flexible, but may require additional customization to work for your viz. Editing them is fine, but note that `'src/components/ReferencesSection.vue'` and `'src/components/AuthorshipSection.vue'` make use of the `'VizSection.vue'` template, so please do not make breaking changes. If you think your edits could be useful in other sites, please submit a MR to the `vue3-template` repo.
6. Development conventions/best practices
    * Upload data files (e.g., `.csv` or `.json` files) svgs and images to s3 - don't commit them to the repository.
    * Page styling with `css`
      * Import all fonts and set **global css color variables** and body and html element styling (including **global font sizing**) in `'src/assets/css/base.css'`
        * Note that we use a [fixed `:root` font size of `62.5%`](https://blog.hubspot.com/website/css-rem#:~:text=By%20setting%20the%20root%20font,%2C%20and%202.0rem%2C%20respectively.), and all `css` font sizes should be specified in units of `rem`, where `1 rem` is equivalent `10 px`.
      * Use `'main.css'` for global page content styling that is exclusive of the USGS header and footer and is not component-specific (e.g., defining styles for standard text or figure container `<div>` elements, styling all section titles, etc.). For all colors, remember to reference color variables set in `'src/assets/css/base.css'`, e.g., `color: var(--color-title-text);`
      * Put component-specific styling in the `<style>` tags of specific components. Again, for all colors, remember to reference color variables set in `'src/assets/css/base.css'`, e.g., `color: var(--color-title-text);`
    * Class and ID naming conventions
      * TBD
      * Do not use spaces in class or ID names
    * JavaScript conventions
      * TBD  
