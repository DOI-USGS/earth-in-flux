<template>
    <section>
        <!---VizSection-->
        <VizSection
            :figures="false"
            :fig-caption="false"
        >
            <!-- HEADING -->
            <template #heading>
                <h2 v-html="text.heading1" />
            </template>
            <!-- FIGURES -->
            <template #aboveExplanation>
                <p v-html="text.paragraph1" />
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
                <p v-html="tab.tabText" v-if="primaryCategorySelected"/>
                <p v-html="subThreatText" v-if="!primaryCategorySelected"/>
                <div id="icon-legend-container">
                    <img class="tab-icon-image" :src="iconSource" alt="">
                    <img class="tab-legend-image" :src="legendSource" :alt="tab.tabLegendImageAlt">
                </div>
                <img class="tab-map-image" :src="mapSource" :alt="tab.tabMapImageAlt">
            </tabItem>
        </tabsGroup>
    </section>
</template>

<script setup>
    import { computed, nextTick, onMounted, ref } from "vue";
    import VizSection from '@/components/VizSection.vue';
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
    const subThreatText = computed(() => {
        return subCategoryData.value.subThreatText;
    });

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

</script>

<style lang="scss">
$habitat: #4E6D6E; 
$habitat-faded: #C9D8D9;
$habitat-dark: #405959;
$pollution: #7A562B;
$pollution-faded: #E1C8AA;
$pollution-dark: #5B401F;
$climate-and-weather: #002D5E;
$climate-and-weather-faded: #B2C0CE;
$climate-and-weather-dark: #002D5E;
$invasive-species: #B74F49;
$fishing-pressure: #835192;
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
.tab-map-image {
    width: 100%;
}
.habitat {
    color: $habitat;
    fill: $habitat;
}
.pollution {
    color: $pollution;
    fill: $pollution;
}
.climate-and-weather {
    color: $climate-and-weather;
    fill: $climate-and-weather;
}
.invasive-species {
    color: $invasive-species;
    fill: $invasive-species;
}
.fishing-pressure {
    color: $fishing-pressure;
    fill: $fishing-pressure;
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
.highlight.fishing-pressure {
    background-color: $fishing-pressure;
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
        background-color: $habitat-faded;
        color: $habitat-dark;
    }
    .category-button.pollution:hover {
        background-color: $pollution-faded;
        color: $pollution-dark;
    }
    .category-button.climate-and-weather:hover {
        background-color: $climate-and-weather-faded;
        color: $climate-and-weather-dark;
    }
    .category-button.invasive-species:hover {
        color: white;
    }
    .category-button.fishing-pressure:hover {
        color: white;
    }
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