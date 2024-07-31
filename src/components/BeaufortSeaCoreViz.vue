<template>
    <!---VizSection-->
    <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div class="single maxWidth">
                <figure>
                    <img src="../../beaufortSea/images/BeaufortSeaCore_1.png">
                </figure>
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            <p v-html="text.paragraph1" />
        </template>
    </VizSection>
        <!---VizSection-->
        <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div class="single maxWidth">
                <figure>
                    <img src="../../beaufortSea/images/BeaufortSeaCore_2.png">
                </figure>
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            <p v-html="text.paragraph2" />
        </template>
    </VizSection>
        <!---VizSection-->
        <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div class="single maxWidth">
                <figure>
                    <img src="../../beaufortSea/images/BeaufortSeaCore_3.gif">
                </figure>
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            <p v-html="text.paragraph3" />
        </template>
    </VizSection>
        <!---VizSection-->
        <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div class="single maxWidth">
                <figure>
                    <img src="../../beaufortSea/images/BeaufortSeaCore_4.png">
                </figure>
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            <p v-html="text.paragraph4" />
        </template>
    </VizSection>
        <!---VizSection-->
        <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div class="single maxWidth">
                <figure>
                    <img src="../../beaufortSea/images/BeaufortSeaCore_5.gif">
                </figure>
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            <p v-html="text.paragraph5" />
        </template>
    </VizSection>
        <!---VizSection-->
        <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div class="single maxWidth">
                <figure>
                    <img src="../../beaufortSea/images/BeaufortSeaCore_6.png">
                </figure>
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            <p v-html="text.paragraph6" />
        </template>
    </VizSection>
        <!---VizSection-->
        <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div class="single maxWidth">
                <figure>
                    <img src="../../beaufortSea/images/BeaufortSeaCore_7.gif">
                </figure>
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            <p v-html="text.paragraph7" />
        </template>
    </VizSection>
        <!---VizSection-->
        <VizSection
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <!-- FIGURES -->
        <template #aboveExplanation>
        </template>
        <template #figures>
            <div class="single maxWidth">
                <figure>
                    <img src="../../beaufortSea/images/BeaufortSeaCore_8.png">
                </figure>
            </div>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
            <p v-html="text.paragraph8" />
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