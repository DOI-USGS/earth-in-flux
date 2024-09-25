<template>
    <div>
        <!---VizSection-->
        <VizSection
            :figures="false"
            :fig-caption="false"
        >
            <template #heading>
                <h2>
                    {{ text.heading1 }}
                </h2>
            </template>
            <template #aboveExplanation>
                <p v-html="text.introParagraph1" />
                <p v-html="text.introParagraph2" />
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
        <!---VizSection 2-->
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <!-- HEADING -->
            <template #heading>
                <h2>
                    {{ text.heading2 }}
                </h2>
            </template>
            <!-- FIGURES -->
            <template #aboveExplanation>
            </template>
            <template #figures>
                <div class="single narrow">
                    <figure>
                        <img src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaSpecies_2.png">
                    </figure>
                </div>
            </template>
            <!-- EXPLANATION -->
            <template #belowExplanation>
                {{ text.paragraph2 }}
            </template>
        </VizSection>
        <!---VizSection2a-->
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <!-- HEADING -->
            <template #heading>
                <h3 v-html="text.heading2a" />
            </template>
            <!-- FIGURES -->
            <template #aboveExplanation>
                <p v-html="text.paragraph2a" />
            </template>
            <template #figures>
                <div class="single narrow">
                    <figure>
                        <img src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaSpecies_2a.png">
                    </figure>
                </div>
            </template>
            <!-- EXPLANATION -->
            <template #belowExplanation>
            </template>
        </VizSection>
        <!---VizSection2b-->
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <template #heading>
                <h3 v-html="text.heading2b" />
            </template>
            <!-- FIGURES -->
            <template #aboveExplanation>
                <p v-html="text.paragraph2b" />
            </template>
            <template #figures>
                <div class="single narrow">
                    <figure>
                        <img src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaSpecies_2b.png">
                    </figure>
                </div>
            </template>
            <!-- EXPLANATION -->
            <template #belowExplanation>
            </template>
        </VizSection>
        <!---VizSection 2c-->
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <!-- HEADING -->
            <template #heading>
                <h3 v-html="text.heading2c" />
            </template>
            <!-- FIGURES -->
            <template #aboveExplanation>
                <p v-html="text.paragraph2c" />
            </template>
            <template #figures>
                <div class="single narrow">
                    <figure>
                        <img src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaSpecies_2c.png">
                    </figure>
                </div>
            </template>
            <!-- EXPLANATION -->
            <template #belowExplanation>
            </template>
        </VizSection>
        <!---VizSection 3-->
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <!-- HEADING -->
            <template #heading>
                <h2>
                    {{ text.heading3 }}
                </h2>
            </template>
            <!-- FIGURES -->
            <template #aboveExplanation>
            </template>
            <template #figures>
                <div class="single narrow">
                    <figure>
                        <img src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaSpecies_3.png">
                    </figure>
                </div>
            </template>
            <!-- EXPLANATION -->
            <template #belowExplanation>
                <p v-html="text.paragraph3" />
            </template>
        </VizSection>
        <!---VizSection 3a-->
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <!-- HEADING -->
            <template #heading>
                <h3 v-html="text.heading3a" />
            </template>
            <!-- FIGURES -->
            <template #aboveExplanation>
                <p v-html="text.paragraph3a" />
            </template>
            <template #figures>
                <div class="single narrow">
                    <figure>
                        <img src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaSpecies_3a.png">
                    </figure>
                </div>
            </template>
            <!-- EXPLANATION -->
            <template #belowExplanation>
            </template>
        </VizSection>
        <!---VizSection 3b-->
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <!-- HEADING -->
            <template #heading>
                <h3 v-html="text.heading3b" />
            </template>
            <!-- FIGURES -->
            <template #aboveExplanation>
                <p v-html="text.paragraph3b" />
            </template>
            <template #figures>
                <div class="single narrow">
                    <figure>
                        <img src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaSpecies_3b.png">
                    </figure>
                </div>
            </template>
            <!-- EXPLANATION -->
            <template #belowExplanation>
            </template>
        </VizSection>
    </div>
</template>

<script setup>
    // import {Tabs, Tab} from 'vue3-tabs-component';
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })

    function getContentImageUrl(suffix) {
        return `https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaSpecies_${suffix}.png`
    }
    function getPrefixImageHTML(image_name) {
        return `<img class='tab-image' src='https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/${image_name}.png'>`
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
.tabs-component {
    margin: auto;
    width: 90vw;
}

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
.tab-image {
    max-width: 5rem;
    max-height: 5rem;
    margin-right: 1rem;
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
</style>