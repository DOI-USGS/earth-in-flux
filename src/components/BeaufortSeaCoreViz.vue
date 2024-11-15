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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_1.png" alt="">
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_2.png" alt="">
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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_3.gif" alt="">
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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_4.png" alt="">
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_5.gif" alt="">
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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_6.png" alt="">
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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_7.gif" alt="">
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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_8.png" alt="">
                </div>
            </template>
        </VizSection>
    </div>
</template>

<script setup>
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })
</script>

<style>
.coring-image-container {
    margin: 3rem auto 4rem auto;
}
.image-narrow {
    max-width: 80rem;
    justify-self: center;
}    
.coring-image {
        width: 100%;
        justify-self: center;
        @media only screen and (max-width: 600px) {
            width: 100%;
        }
    }
</style>