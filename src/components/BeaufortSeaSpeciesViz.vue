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
.subheading-container {
    margin: 1rem 0 1rem 0;
    height: 5rem;
}
.subheading-image {
    height: 5rem;
    margin: 0 1rem 0 1rem;
}
.subheading {
    padding: 0;
    display: inline-block;
    transform: translate(0%, -50%);
}
.tab-content-image {
    width: 100%;
}
.tab-content-title {
    padding: 0rem 0 1rem 0;
    line-height: 2.6rem;
}
</style>