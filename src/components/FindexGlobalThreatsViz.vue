<template>
    <section>
        <!---VizSection-->
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <template #heading>
                <h2 v-html="text.heading1" />
            </template>
            <template #aboveExplanation>
                <p v-html="text.paragraph1" />
            <div class="toggle-group">
                <p id="toggle-title">Toggle map layers:</p> 
                <ToggleSwitch 
                    v-for="layer in layers"
                    :key="layer.order"
                    v-model="layer.visible" 
                    :label="layer.label"
                    :rightColor="layer.color"
                />
            </div>
            </template>
            <template #figures>    
                <div id="top-threat-map-container">
                    <figure
                        class="base-image"
                    >
                        <img class="map-image" :src="getImageURL('base_threat_by_basin.png')" alt="">
                    </figure>
                    <figure 
                        v-for="layer in layers"
                        :key="layer.order"
                        v-show="layer.visible"
                        class="overlay-image"
                    >
                        <img class="map-image" :src="getImageURL(layer.path)" alt="">
                    </figure>
                </div>
            </template>
        </VizSection>
        <!---VizSection-->
        <VizSection
            :figures="false"
            :fig-caption="false"
        >
            <template #heading>
                <h2 v-html="text.heading2" />
            </template>
            <template #aboveExplanation>
                <p v-html="text.paragraph2" />
            </template>
        </VizSection>
        <tabsGroup id="map-tabs" :options="{ useUrlFragment: false }" >
            <tabItem v-for="tab in text.tabData" :name="tab.tabTitle" :key="tab.tabTitle" :prefix="getPrefixImageHTML(tab.tabIcon)">
                <h3>
                    Threat category:
                    <button v-html="tab.tabContentTitle" @click = "switchToPrimaryCategory(tab.tabContentTitle, tab.subThreatPrefix)" :class="[tab.tabContentTitleID, { 'highlight': currentCategory == tab.tabContentTitle }]" :id="tab.tabContentTitleID" class="category-button primary" />
                </h3>
                <div id="button-container">
                    <h4 v-if="tab.subThreatData.length > 1">Subthreat categories:
                        <span v-for="subThreatCategory, index in tab.subThreatData" :key="subThreatCategory.subThreat">
                            <button @click="switchToSubCategory(subThreatCategory.subThreat, tab.subThreatPrefix)" :class="[tab.tabContentTitleID, { 'highlight': currentCategory == subThreatCategory.subThreat }]" v-html="subThreatCategory.subThreat" class="category-button sub"></button>
                            <span v-if="index < tab.subThreatData.length-1" :class="tab.tabContentTitleID" class="separator">
                                <FishIcon id="findex-fish" :class="tab.tabContentTitleID"/>
                            </span>
                        </span>
                    </h4>
                </div>
                <div id="icon-legend-container">
                    <img class="tab-icon-image" :src="iconSource" alt="">
                    <img class="tab-legend-image" :src="legendSource" :alt="tab.tabLegendImageAlt">
                </div>
                <div id="threat-map-container">
                    <img class="tab-map-image" :src="mapSource" :alt="tab.tabMapImageAlt">
                </div>
                <p v-html="tab.tabText" v-if="primaryCategorySelected"/>
                <div v-if="!primaryCategorySelected">
                    <CollapsibleAccordion 
                        v-for="item, index in subCategoryData.subThreatText"
                        :key="index"
                        :heading="item.heading"
                        :content="item.content"
                        :active-on-load="item.activeOnLoad"
                        :left-border-color="tab.tabColor"
                        button-active-background-color="var(--light-grey)"
                        button-inactive-background-color="var(--light-grey)"
                        button-font-weight="bold"
                        button-font-color="var(--color-text)"
                    />
                </div>
            </tabItem>
        </tabsGroup>
    </section>
</template>

<script setup>
    import { computed, nextTick, onMounted, reactive, ref } from "vue";
    import VizSection from '@/components/VizSection.vue';
    import ToggleSwitch from '@/components/ToggleSwitch.vue';
    import CollapsibleAccordion from '@/components/CollapsibleAccordion.vue';
    import FishIcon from '@/assets/svgs/noun-fish-7471722.svg';

    // define props
    const props = defineProps({
        text: { type: Object }
    })

    // Global variables 
    let primaryCategorySelected = ref(true);
    let currentCategory = ref(null);
    let currentCategorySubThreatPrefix = ref(null);
    let legendSource = ref(null);
    let mapSource = ref(null);
    let iconSource = ref(null);

    // Set up computed variables that depend on ref values
    const primaryCategoryData = computed(() => {
        return props.text.tabData.filter(d => d.subThreatPrefix == currentCategorySubThreatPrefix.value)[0]
    })
    // undefined if primaryCategorySelected
    const subCategoryData = computed(() => {
        return primaryCategoryData.value.subThreatData.filter(d => d.subThreat == currentCategory.value)[0]
    })

    // set up reactive layers object
    const layers = reactive(props.text.mapData)

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        try {            
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
        const activeTab = tabs.querySelectorAll(".is-active a")
        
        // pull category information
        currentCategory.value = activeTab[0].text
        const currentData = props.text.tabData.filter(d => d.tabContentTitle == currentCategory.value)[0]
        currentCategorySubThreatPrefix.value = currentData.subThreatPrefix

        // update map - always show primary category on page load
        switchToPrimaryCategory(currentCategory.value, currentCategorySubThreatPrefix.value)
    }
    function getPrefixImageURL(filename) {
        return new URL(`../assets/svgs/${filename}.svg`, import.meta.url).href
    }
    function getPrefixImageHTML(filename) {
        const imgURL = getPrefixImageURL(filename)
        return `<img class='tab-image' src=${imgURL}>`
    }
    function getContentImageUrl(title, category_prefix, content_type) {
        if (primaryCategorySelected.value) {
            return new URL(`../assets/images/${title.replace(/ /g, "_")}_${content_type}.png`, import.meta.url).href
        } else {
            return new URL(`../assets/images/${category_prefix}_${title.replace(/ /g, "_")}_${content_type}.png`, import.meta.url).href
        }        
    }
    function updateTabContent(category, prefix) {
        currentCategory.value = category
        currentCategorySubThreatPrefix.value = prefix
        mapSource.value = getContentImageUrl(currentCategory.value, currentCategorySubThreatPrefix.value, "map")
        legendSource.value = getContentImageUrl(currentCategory.value, currentCategorySubThreatPrefix.value, "legend")
    }
    function updateIcon() {
        if (primaryCategorySelected.value) {
            iconSource.value = getPrefixImageURL(primaryCategoryData.value.tabIcon)
        } else {
            iconSource.value = getPrefixImageURL(subCategoryData.value.subThreatIcon)
        } 
    }
    function switchToSubCategory(category, prefix) {
        primaryCategorySelected.value = false;
        updateTabContent(category, prefix);
        updateIcon();
    }
    function switchToPrimaryCategory(category, prefix) {
        primaryCategorySelected.value = true;
        updateTabContent(category, prefix);
        updateIcon();
    }
    function getImageURL(file) {
        return new URL(`../assets/images/${file}`, import.meta.url).href
    }

</script>

<style lang="scss">
#toggle-title {
    font-weight: 700;
}
#top-threat-map-container {
    position: relative;
    margin: 4rem auto 3rem auto;
    max-width: 1200px;
}
.base-image {
    position: relative;
}
.overlay-image {
    position: absolute;
    top: 0;
    left: 0;
}
.map-image {
    width: 100%;
    justify-self: center;
    @media only screen and (max-width: 600px) {
        width: 100%;
    }
}
#map-tabs {
    margin-top: 3rem;
}
#button-container {
    margin-bottom: 2rem;
}
#findex-fish {
    width: 20px;
    height: 20px;
    transform: translate(0, 4px);
}
#icon-legend-container{
    display: flex;
    flex-direction: row;
    align-items: center;
}
.tab-icon-image {
    max-width: 50px;
    max-height: 50px;
    margin-right: 1rem;
    height: auto;
    width: auto;
    @media only screen and (max-width: 600px) {
        max-width: 40px;
        max-height: 40px;
    }
}
.tab-legend-image {
    width: 230px;
    @media only screen and (max-width: 600px) {
        width: min(190px, 100vw);
    }
}
#threat-map-container {
    margin: 4rem 0 5rem 0;
}
.tab-map-image {
    width: 100%;
}
.habitat {
    color: var(--color-habitat);
    fill: var(--color-habitat);
}
.pollution {
    color: var(--color-pollution);
    fill: var(--color-pollution);
}
.climate-and-weather {
    color: var(--color-climate-and-weather);
    fill: var(--color-climate-and-weather);
}
.invasive-species {
    color: var(--color-invasive-species);
    fill: var(--color-invasive-species);
}
.fishing-pressure {
    color: var(--color-fishing-pressure);
    fill: var(--color-fishing-pressure);
}
.highlight.habitat {
    background-color: var(--color-habitat);
}
.highlight.pollution {
    background-color: var(--color-pollution);
}
.highlight.climate-and-weather {
    background-color: var(--color-climate-and-weather);
}
.highlight.invasive-species {
    background-color: var(--color-invasive-species);
}
.highlight.fishing-pressure {
    background-color: var(--color-fishing-pressure);
}
.category-button {
    background-color: transparent;
    border: 0rem;
    padding: 0.05rem 0.8rem 0.2rem 0.75rem;
    border-radius: 10px;
    text-decoration: underline;
}
.category-button.highlight {
    text-decoration: none;
    color: white;
}
@media (hover: hover) {
    .category-button.habitat:hover {
        background-color: var(--color-habitat-faded);
        color: var(--color-habitat-dark);
    }
    .category-button.pollution:hover {
        background-color: var(--color-pollution-faded);
        color: var(--color-pollution-dark);
    }
    .category-button.climate-and-weather:hover {
        background-color: var(--color-climate-and-weather-faded);
        color: var(--color-climate-and-weather-dark);
    }
    .category-button.invasive-species:hover {
        color: white;
    }
    .category-button.fishing-pressure:hover {
        color: white;
    }
}
.habitat {
    text-decoration: underline solid var(--color-habitat-faded);
}
.pollution {
    text-decoration: underline solid var(--color-pollution-faded);
} 
.climate-and-weather {
    text-decoration: underline solid var(--color-climate-and-weather-faded);
} 
.separator {
    background-color: transparent;
    text-decoration: none;
}
</style>