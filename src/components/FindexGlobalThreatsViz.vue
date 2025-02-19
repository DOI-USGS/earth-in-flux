<template>
    <!---VizSection-->
    <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <template #heading>
            <h2>
            </h2>
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
                <p v-html="text.paragraph1" />
        </template>
        <template #figures>
            <tabsGroup id="map-tabs" :options="{ useUrlFragment: false }" >
                <tabItem v-for="tab in text.tabData" :name="tab.tabTitle" :key="tab.tabTitle" :prefix="getPrefixImageHTML(tab.tabIcon)">
                    <h3>
                        Threat category:
                        <button v-html="tab.tabContentTitle" @click = "switchToPrimaryCategory(tab.tabContentTitle, tab.subThreatPrefix)" :class="[tab.tabContentTitleID, { 'active': currentCategory == tab.tabContentTitle, 'highlight': currentCategory == tab.tabContentTitle }]" :id="tab.tabContentTitleID" class="category-button" />
                    </h3>
                    <div id="button-container">
                        <h4 v-if="tab.subThreatData.length > 1">Subthreat categories:
                            <span v-for="subThreat, index in tab.subThreatData" :key="subThreat">
                                <button @click="switchToSubCategory(subThreat, tab.subThreatPrefix)" :class="[tab.tabContentTitleID, { 'active': currentCategory == subThreat, 'highlight': currentCategory == subThreat }]" v-html="subThreat" class="subcategory-button"></button>
                                <span v-if="index < tab.subThreatData.length-1"> | </span>
                            </span>
                        </h4>
                    </div>
                    <p v-html="tab.tabText" />
                    <img class="tab-content-image" :src="mapSource" :alt="tab.tabImageAlt">
                </tabItem>
            </tabsGroup>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
        </template>
    </VizSection>
</template>

<script setup>
    import { nextTick, onMounted, ref } from "vue";
    // import { isMobile } from 'mobile-device-detect';
    import VizSection from '@/components/VizSection.vue';

    // define props
    const props = defineProps({
        text: { type: Object }
    })

    // Global variables 
    // const publicPath = import.meta.env.BASE_URL;
    // const mobileView = isMobile;
    let primaryCategorySelected = ref(true);
    let currentCategory = ref(null);
    let currentCategorySubThreatPrefix = ref(null);
    let mapSource = ref(null);

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        try {
            console.log("on mounted")
            
            // Once DOM is up to date, make sure the currently shown tab is updated
            await nextTick()
            updateTab()

            // add event listeners to all tabs, so that update on click
            const tabs = document.getElementById("map-tabs")
            const tabList = tabs.querySelectorAll("a")
            tabList.forEach(tab => {
                tab.href = ""
                tab.addEventListener("click", updateTab)
            })
        } catch (error) {
            console.error('Error during component mounting', error);
        }        
    });
    function updateTab() {
        // identify active tab
        const tabs = document.getElementById("map-tabs")
        console.log(tabs)
        console.log(tabs.querySelectorAll("a"))
        const activeTab = tabs.querySelectorAll(".is-active a")
        console.log(activeTab)
        
        // pull category information
        currentCategory.value = activeTab[0].text
        const currentData = props.text.tabData.filter(d => d.tabContentTitle == currentCategory.value)
        currentCategorySubThreatPrefix.value = currentData.subThreatPrefix

        // update map
        switchToPrimaryCategory(currentCategory.value, currentCategorySubThreatPrefix.value)
    }
    function getPrefixImageURL(filename) {
        return `src/assets/images/${filename}.png`
    }
    function getPrefixImageHTML(filename) {
        const imgURL = getPrefixImageURL(filename)
        return `<img class='tab-image' src=${imgURL}>`
    }
    function getContentImageUrl(title, category_prefix) {
        if (primaryCategorySelected.value) {
            return `src/assets/images/${title.replace(/ /g, "_")}_map.png`
        } else {
            return `src/assets/images/${category_prefix}${title.replace(/ /g, "_")}_map.png`
        }        
    }
    function switchToSubCategory(category, prefix) {
        primaryCategorySelected.value = false;
        currentCategory.value = category
        currentCategorySubThreatPrefix.value = prefix
        console.log(`switching to ${currentCategory.value}, ${currentCategorySubThreatPrefix.value}`)
        mapSource.value = getContentImageUrl(currentCategory.value, currentCategorySubThreatPrefix.value)
    }
    function switchToPrimaryCategory(category, prefix) {
        primaryCategorySelected.value = true;
        currentCategory.value = category
        currentCategorySubThreatPrefix.value = prefix
        mapSource.value = getContentImageUrl(currentCategory.value, currentCategorySubThreatPrefix.value)
    }

</script>

<style>

#map-tabs {
    margin-top: 3rem;
}
#button-container {
    margin-bottom: 2rem;
}
.active {
        /* border-color: ; */
        border-width: 3px;
        border-style: solid;
        /* border-radius: 1.5px; */
        opacity: 1;
        font-weight: 500;
        color: white;
        /* transform: scale(1.1); */
        
    }
/* .active.habitat {
    background-color: #7A562B;
} */
.category-button {
    background-color: transparent;
    padding-left: 0.75rem;
    padding-bottom: 0.2rem;
    border: 0rem;
    padding: 0.05rem 0.8rem 0.2rem 0.75rem;
    border-radius: 10px;
    white-space: nowrap;
    font-weight: bold;
    font-size: 2.2rem;
    margin-bottom: 1rem;
}
.category-button:hover {
    color: var(--color-text)
}
.category-button:hover.habitat {
    background-color: #E1C8AA;
}
.category-button:hover.pollution {
    background-color: #B2C0CE;
}
.category-button:hover.climate-and-weather {
    background-color: #DDCCE2;
}
.category-button:hover.invasive-species {
    color: white;
    /* background-color: #C9D8D9; */
}
.category-button:hover.exploitation {
    color: white;
    /* background-color: #E1C8AA; */
}
.highlight.habitat {
    background-color: #7A562B;
}
.highlight.pollution {
    background-color: #002D5E;
}
.highlight.climate-and-weather {
    background-color: #835192;
}
.highlight.invasive-species {
    background-color: #4E6D6E;
}
.highlight.exploitation {
    background-color: #B74F49;
}
.subcategory-button {
    background-color: transparent;
    border: 0rem;
    border-radius: 10px;
    padding: 0.05rem 0.8rem 0.2rem 0.75rem;
}
.subcategory-button:hover {
    color: var(--color-text)
}
.subcategory-button:hover.habitat {
    background-color: #E1C8AA;
}
.subcategory-button:hover.pollution {
    background-color: #B2C0CE;
}
.subcategory-button:hover.climate-and-weather {
    background-color: #DDCCE2;
}
ul {
    padding-inline-start: 0px;
}
li {
    padding: 0; 
}
.tabs-component {
    margin: auto;
    width: 90vw;
    max-width: 1000px;
}
@media (min-width: 1000px) {
    .tabs-component {
        width: 70vw;
    }
}
.tab-image {
    max-width: 2.5rem;
    max-height: 2.5rem;
    margin-right: 1rem;
    height: auto;
    width: auto;
}

@media (min-width: 1000px) {
    .tab-image {
        max-width: fit-content;
        max-height: 5rem;
        margin-bottom: 1rem;
        height: 5rem;
        width: auto;
    }
}
.subheading-container {
    margin: 1rem 0 1rem 0;
    height: 5rem;
}
.subheading-image {
    /* max-width: 5rem; */
    height: 5rem;
    margin: 0 1rem 0 1rem;
}
.subheading {
    padding: 0;
    display: inline-block;
    transform: translate(0%, -50%);
}
.tabs-component-tab.is-inactive img {
    opacity: 0.5;
}
.tab-content-image {
    width: 100%;
}
.tab-content-title {
    padding: 1rem 1rem 1rem 1rem;
    line-height: 2.6rem;
}
/* css modified from vue-tabs-component demo: https://tabs-component.jakubpotocky.sk/ */
.tabs-component-tabs {
    border: solid 1px #ddd;
    border-radius: 6px;
    margin-bottom: 5px;
}

@media (min-width: 1000px) {
    .tabs-component-tabs {
        border: 0;
        align-items: stretch;
        display: flex;
        justify-content: flex-start;
        margin-bottom: -1px;
    }
}

.tabs-component-tab, .tabs-component-tab--custom {
    color: #999;
    font-size: 1.6rem;
    font-weight: 600;
    margin-right: 0;
    list-style: none;
}

.tabs-component-tab:not(:last-child) {
    border-bottom: dotted 1px #ddd;
}

.tabs-component-tab:hover {
    color: #666;
}

.tabs-component-tab.is-active {
    color: #000;
}
.tabs-component-tab.is-disabled *, .tabs-component-tab--custom.is-disabled * {
    color: #cdcdcd;
    cursor: not-allowed !important;
}

@media (min-width: 1000px) {
    .tabs-component-tab, .tabs-component-tab--custom {
        background-color: #fff;
        border: solid 1px #ddd;
        border-radius: 3px 3px 0 0;
        margin-right: .5em;
        /* transform: translateY(4px); */
        transition: transform .3s ease;
        font-size: 1.6rem;
        text-align: center;
    }

    .tabs-component-tab.is-active, .tabs-component-tab--custom.is-active {
        border-bottom: solid 1px #fff;
        z-index: 2;
        transform: translateY(0);
    }
}

@media (min-width: 1000px) {
    .tabs-component-tab-a, .tabs-component-tab-a--custom {
        align-items: center;
        color: inherit;
        display: flex;
        flex-direction: column;
        padding: .75em 1em;
        text-decoration: none;
    }
}

.tabs-component-tab-a, .tabs-component-tab-a--custom {
    align-items: center;
    color: inherit;
    display: flex;
    padding: .75em 1em;
    text-decoration: none;
}

.tabs-component-panels {
    padding: 2em 1em;
    background-color: #fff;
    border: solid 1px #ddd;
    border-radius: 0 6px 6px 6px;
    box-shadow: 0 0 10px rgba(0, 0, 0, .05);
}

@media (min-width: 1000px) {
    .tabs-component-panels {
        background-color: #fff;
        border: solid 1px #ddd;
        border-radius: 0 6px 6px 6px;
        box-shadow: 0 0 10px rgba(0, 0, 0, .05);
        padding: 2em 2em;
    }
}

.tabs-component-btn {
  cursor: pointer;
  background: #e1ecf4;
  border-radius: 3px;
  border: 1px solid #7aa7c7;
  padding: 4px 8px;
  color: #39739d;
}

.tabs-component-btn:hover {
  background-color: #b3d3ea;
  color: #2c5777;
}

.tabs-component-btn:active {
  background-color: #a0c7e4;
  box-shadow: none;
  color: #2c5777;
}

.tabs-component-tab--custom {
    border-color: #e1ecf4;
    color: #68838d;
}

.tabs-component-tab--custom:hover {
    border-color: #e1ecf4;
    color: #39739d;
}

.tabs-component-tab--custom.is-active {
    color: #39739d;
    border-color: #7aa7c7;
    border-bottom: solid 1px #fff;
}
</style>