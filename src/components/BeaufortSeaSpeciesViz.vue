<template>
    <div>
        <!---VizSection-->
        <VizSection
            :figures="false"
            :fig-caption="false"
        >
            <template #heading>
                <h2 v-html="text.heading1" />
            </template>
            <template #aboveExplanation>
                <p v-html="text.paragraph1" />
            </template>
        </VizSection>
        
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <template #heading>
                <h2 v-html="text.heading2" />
            </template>
            <template #aboveExplanation>
                <p v-html="text.paragraph2" />
                <div class="subheading-container">
                    <h3 class="subheading" v-html="text.subheading1" />
                    <img v-for="foram in forams" class='subheading-image' :src="getPrefixImageURL(foram.tabPrefixImageName)" :key="foram.tabTitle"/>
                </div>
                <p v-html="text.paragraph3" />
                <div class="subheading-container">
                    <h3 class="subheading" v-html="text.subheading2" />
                    <img v-for="ostracode in ostracodes" class='subheading-image' :src="getPrefixImageURL(ostracode.tabPrefixImageName)" :key="ostracode.tabTitle"/>
                </div>
                <p v-html="text.paragraph4" />
            </template>
        </VizSection>
        <tabsGroup id="species-tabs" :options="{ useUrlFragment: false }">
            <tabItem v-for="tab in text.tabData" :name="`<span class='scientificName'>${tab.tabTitle}</span>`" :key="tab.tabTitle" :prefix="getPrefixImageHTML(tab.tabPrefixImageName)">
                <h3 class="tab-content-title">
                    <span v-if="!mobileView" class="species-class" :id="`species-${tab.tabContentTitleID}`">
                        <span class="scientificName highlight species-title" :id="tab.tabContentTitleID">
                            {{ tab.tabContentTitle }}
                        </span>
                        {{ tab.tabSpeciesClass }}
                    </span>
                    <span v-if="mobileView" class="scientificName highlight species-title" :id="tab.tabContentTitleID">
                        {{ tab.tabContentTitle }}
                    </span>
                </h3>
                <h3 v-if="mobileView" class="tab-content-title">
                    <span class="species-class mobile" :id="`species-${tab.tabContentTitleID}`">
                        {{ tab.tabSpeciesClass }}
                    </span>
                </h3>
                <p v-html="tab.tabText" />
                <img class="tab-content-image" :src="getContentImageUrl(tab.tabContentImageSuffix)" :alt="tab.tabImageAlt">
            </tabItem>
        </tabsGroup>
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <template #heading>
                <h2 v-html="text.heading3" />
            </template>
            <template #aboveExplanation>
                <p v-html="text.paragraph5" />
                <p v-html="text.paragraph6" />
            </template>
        </VizSection>
    </div>
</template>

<script setup>
    import { isMobile } from 'mobile-device-detect';
    import VizSection from '@/components/VizSection.vue';

    // define props
    const props = defineProps({
        text: { type: Object }
    })

    // global objects
    const mobileView = isMobile;
    const baseURL = "https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/"
    const forams = props.text.tabData.filter(item => item.tabPrefixImageName.includes('F'))
    const ostracodes = props.text.tabData.filter(item => item.tabPrefixImageName.includes('O'))

    function getContentImageUrl(suffix) {
        return baseURL + `BeaufortSeaSpecies_${suffix}.webp`
    }
    function getPrefixImageURL(image_name) {
        return baseURL + `${image_name}.webp`
    }
    function getPrefixImageHTML(image_name) {
        const imgURL = getPrefixImageURL(image_name)
        return `<img class='tab-image' src=${imgURL}>`
    }
</script>

<style>

.species-title {
    padding: 1px 10px 2px 6px;
}
.species-class {
    border-radius: 10px;
    padding: 1px 10px 2px 0px;
    font-weight: 300;
    font-style: italic;
}
.species-class.mobile {
    padding: 1px 10px 2px 6px;
}
#species-cassidulina {
    border: 1px solid #3c475a;
}
#species-elphidium {
    border: 1px solid #66768F;
}
#species-paracyprideis {
    border: 1px solid #729C9D;
}
#species-kotorachythere {
    border: 1px solid #c49051;
}
#species-spiroplectimmina {
    border: 1px solid #dd605a;
}
#species-tabs {
    margin-top: 3rem;
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
    padding: 0rem 0 1rem 0;
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