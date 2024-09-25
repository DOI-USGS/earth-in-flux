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
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })
</script>

<style>
.narrow{
    width: 100%;
    max-width: 720px;
    margin: 0 auto 0 auto;
}
</style>