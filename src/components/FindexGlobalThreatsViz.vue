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
                        <button v-html="tab.tabContentTitle" @click = "switchToPrimaryCategory(tab.tabContentTitle, tab.subThreatPrefix)" :class="[tab.tabContentTitleID, { 'highlight': currentCategory == tab.tabContentTitle }]" :id="tab.tabContentTitleID" class="category-button primary" />
                    </h3>
                    <div id="button-container">
                        <h4 v-if="tab.subThreatData.length > 1">Subthreat categories:
                            <span v-for="subThreat, index in tab.subThreatData" :key="subThreat">
                                <button @click="switchToSubCategory(subThreat, tab.subThreatPrefix)" :class="[tab.tabContentTitleID, { 'highlight': currentCategory == subThreat }]" v-html="subThreat" class="category-button sub"></button>
                                <span v-if="index < tab.subThreatData.length-1" :class="tab.tabContentTitleID" class="separator"> &#x2022; </span>
                            </span>
                        </h4>
                    </div>
                    <p v-html="tab.tabText" />
                    <img class="tab-legend-image" :src="legendSource" :alt="tab.tabLegendImageAlt">
                    <img class="tab-map-image" :src="mapSource" :alt="tab.tabMapImageAlt">
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
    let legendSource = ref(null);
    let mapSource = ref(null);

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
    function getContentImageUrl(title, category_prefix, content_type) {
        if (primaryCategorySelected.value) {
            return `src/assets/images/${title.replace(/ /g, "_")}_${content_type}.png`
        } else {
            return `src/assets/images/${category_prefix}_${title.replace(/ /g, "_")}_${content_type}.png`
        }        
    }
    function switchToSubCategory(category, prefix) {
        primaryCategorySelected.value = false;
        currentCategory.value = category
        currentCategorySubThreatPrefix.value = prefix
        mapSource.value = getContentImageUrl(currentCategory.value, currentCategorySubThreatPrefix.value, "map")
        legendSource.value = getContentImageUrl(currentCategory.value, currentCategorySubThreatPrefix.value, "legend")
    }
    function switchToPrimaryCategory(category, prefix) {
        primaryCategorySelected.value = true;
        currentCategory.value = category
        currentCategorySubThreatPrefix.value = prefix
        mapSource.value = getContentImageUrl(currentCategory.value, currentCategorySubThreatPrefix.value, "map")
        legendSource.value = getContentImageUrl(currentCategory.value, currentCategorySubThreatPrefix.value, "legend")
    }

</script>

<style lang="scss">
$habitat: #7A562B;
$habitat-faded: #E1C8AA;
$habitat-dark: #5B401F;
$pollution: #002D5E;
$pollution-faded: #B2C0CE;
$pollution-dark: #002D5E;
$climate-and-weather: #835192;
$climate-and-weather-faded: #DDCCE2;
$climate-and-weather-dark: #613A69;
$invasive-species: #4E6D6E;
$exploitation: #B74F49;
#map-tabs {
    margin-top: 3rem;
}
#button-container {
    margin-bottom: 2rem;
}
.tab-legend-image {
    width: 230px;
    @media only screen and (max-width: 600px) {
        width: 100%;
    }
}
.tab-map-image {
    width: 100%;
}
.habitat {
    color: $habitat
}
.pollution {
    color: $pollution;
}
.climate-and-weather {
    color: $climate-and-weather;
}
.invasive-species {
    color: $invasive-species;
}
.exploitation {
    color: $exploitation;
}
.highlight.habitat {
    background-color: $habitat;
}
.highlight.pollution {
    background-color: $pollution;
}
.highlight.climate-and-weather {
    background-color: $climate-and-weather;
}
.highlight.invasive-species {
    background-color: $invasive-species;
}
.highlight.exploitation {
    background-color: $exploitation;
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
.category-button:hover {
    color: var(--color-text);
}
.category-button:hover.habitat {
    background-color: $habitat-faded;
    color: $habitat-dark;
}
.category-button:hover.pollution {
    background-color: $pollution-faded;
    color: $pollution-dark;
}
.category-button:hover.climate-and-weather {
    background-color: $climate-and-weather-faded;
    color: $climate-and-weather-dark;
}
.category-button:hover.invasive-species {
    color: white;
}
.category-button:hover.exploitation {
    color: white;
}
.habitat {
    text-decoration: underline solid $habitat-faded;
}
.pollution {
    text-decoration: underline solid $pollution-faded;
} 
.climate-and-weather {
    text-decoration: underline solid $climate-and-weather-faded;
} 
.separator {
    background-color: transparent;
    text-decoration: none;
}
</style>