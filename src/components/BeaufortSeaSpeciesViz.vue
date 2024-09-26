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
                <p v-html="text.paragraph5" />
            </template>
        </VizSection>
        <tabsGroup id="species-tabs" :options="{ useUrlFragment: false }">
            <tabItem v-for="tab, index in text.tabData" :name="tab.tabTitle" :key="tab.tabTitle" :prefix="getPrefixImageHTML(tab.tabPrefixImageName)">
                <h3 class="tab-content-title" v-html="tab.tabTitle" />
                <p v-html="tab.tabText" />
                {{ index }}
                <img class="tab-content-image" :src="getContentImageUrl(tab.tabContentImageSuffix)">
            </tabItem>
        </tabsGroup>
    </div>
</template>

<script setup>
    // import {Tabs, Tab} from 'vue3-tabs-component';
    import VizSection from '@/components/VizSection.vue';

    // define props
    const props = defineProps({
        text: { type: Object }
    })

    // global objects
    const baseURL = "https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/"
    const forams = props.text.tabData.filter(item => item.tabPrefixImageName.includes('F'))
    const ostracodes = props.text.tabData.filter(item => item.tabPrefixImageName.includes('O'))

    function getContentImageUrl(suffix) {
        return baseURL + `BeaufortSeaSpecies_${suffix}.png`
    }
    function getPrefixImageURL(image_name) {
        return baseURL + `${image_name}.png`
    }
    function getPrefixImageHTML(image_name) {
        const imgURL = getPrefixImageURL(image_name)
        return `<img class='tab-image' src=${imgURL}>`
    }
</script>

<style>
.narrow{
    width: 100%;
    max-width: 720px;
    margin: 0 auto 0 auto;
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
    width: 90rem;
}
.tab-image {
    max-width: 5rem;
    max-height: 5rem;
    margin-right: 1rem;
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
}
/* css modified from vue-tabs-component demo: https://tabs-component.jakubpotocky.sk/ */
.tabs-component-tabs {
    border: solid 1px #ddd;
    border-radius: 6px;
    margin-bottom: 5px;
}

@media (min-width: 700px) {
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
    font-size: 14px;
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

@media (min-width: 700px) {
    .tabs-component-tab, .tabs-component-tab--custom {
        background-color: #fff;
        border: solid 1px #ddd;
        border-radius: 3px 3px 0 0;
        margin-right: .5em;
        /* transform: translateY(4px); */
        transition: transform .3s ease;
    }

    .tabs-component-tab.is-active, .tabs-component-tab--custom.is-active {
        border-bottom: solid 1px #fff;
        z-index: 2;
        transform: translateY(0);
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
    padding: 2em 0;
}

@media (min-width: 700px) {
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