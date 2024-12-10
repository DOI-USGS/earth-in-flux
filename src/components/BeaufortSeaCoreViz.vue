<template>
    <div>
        <VizSection
            :figures="false"
            :fig-caption="false"
        >
            <template #heading>
                <h2 v-html="text.heading1" />
            </template>
            <template #aboveExplanation>
                <p v-html="text.intro1" />
                <p v-html="text.intro2" />
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
                <p v-html="text.paragraph1" />
            </template>
            <template #figures>
                <div class="coring-image-container maxWidth group two">
                    <figure v-for="index in getRange(1,2)" :key="index">
                        <img class="coring-image" :src="getImagePath(index)" :alt="getAltText(index)">
                    </figure>
                </div>
            </template>
        </VizSection>
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <template #aboveExplanation>
                <p v-html="text.paragraph2" />
            </template>
            <template #figures>
                <div class="coring-image-container image-narrow">
                    <figure v-for="index in getRange(3,3)" :key="index">
                        <img class="coring-image" :src="getImagePath(index)" :alt="getAltText(index)">
                    </figure>
                </div>
            </template>
        </VizSection>
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <template #aboveExplanation>
                <p v-html="text.paragraph3" />
            </template>
            <template #figures>
                <div class="coring-image-container maxWidth group two">
                    <figure v-for="index in getRange(4,5)" :key="index">
                        <img class="coring-image" :src="getImagePath(index)" :alt="getAltText(index)">
                    </figure>
                </div>
            </template>
        </VizSection>
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <template #aboveExplanation>
                <p v-html="text.paragraph4" />
            </template>
            <template #figures>
                <div class="coring-image-container image-narrow">
                    <figure v-for="index in getRange(6,6)" :key="index">
                        <img class="coring-image" :src="getImagePath(index)" :alt="getAltText(index)">
                    </figure>
                </div>
            </template>
        </VizSection>
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <template #aboveExplanation>
                <p v-html="text.paragraph5" />
            </template>
            <template #figures>
                <div class="coring-image-container image-narrow">
                    <figure v-for="index in getRange(7,7)" :key="index">
                        <img class="coring-image" :src="getImagePath(index)" :alt="getAltText(index)">
                    </figure>
                </div>
            </template>
        </VizSection>
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <template #aboveExplanation>
                <p v-html="text.paragraph6" />
            </template>
            <template #figures>
                <div class="coring-image-container image-narrow">
                    <figure v-for="index in getRange(8,8)" :key="index">
                        <img class="coring-image" :src="getImagePath(index)" :alt="getAltText(index)">
                    </figure>
                </div>
            </template>
        </VizSection>
    </div>
</template>

<script setup>
    import VizSection from '@/components/VizSection.vue';

    // define props
    const props = defineProps({
        text: { type: Object }
    })

    // global variables
    const gifImageIndices = [3, 5, 7]

    function getRange(start, end) {
        return Array.from({ length: end - start + 1 }, (_, i) => i + start)
    }

    function getImagePath(index) {
        let fileEnding;
        if (gifImageIndices.includes(index)) {
            fileEnding = 'gif';
        } else {
            fileEnding = 'webp'
        }
        return `https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_${index}.${fileEnding}`;
    }

    function getAltText(index) {
        return props.text[`alt${index}`]
    }
</script>

<style>
.coring-image-container {
    margin: 3rem auto 4rem auto;
}
.coring-image {
    width: 100%;
    justify-self: center;
    @media only screen and (max-width: 600px) {
        width: 100%;
    }
}
#cesium-tooltip {
    margin-left: -155px;
}
#lead-tooltip {
    margin-left: -145px;
}
#radiocarbon-tooltip {
    margin-left: -185px;
}
</style>