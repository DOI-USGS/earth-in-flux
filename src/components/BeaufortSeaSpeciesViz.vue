<template>
    <!---VizSection-->
    <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <template #heading>
            <h2>
                {{ text.heading1 }}
            </h2>
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div class="single maxWidth">
                <figure>
                    <img src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaSpecies_1.png">
                </figure>
            </div>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            {{ text.paragraph1 }}
        </template>
    </VizSection>
    <!---VizSection-->
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
            <div class="single maxWidth">
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
    <!---VizSection-->
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
            <div class="single maxWidth">
                <figure>
                    <img src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaSpecies_3.png">
                </figure>
            </div>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            {{ text.paragraph3 }}
        </template>
    </VizSection>
</template>

<script setup>
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })
</script>

<style>
</style>