<template>
    <!---VizSection-->
    <VizSection
        id="cross-section"
        :figures="true"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <template #heading>
            <h2>
                {{ text.heading }}
            </h2>
        </template>
        <!-- FIGURES -->
        <template #aboveExplanation>
            <p v-html="text.paragraph1" />
            <p v-html="text.paragraph2" />
            <p v-html="text.paragraph3" />
            <p v-html="text.paragraph4" />
        </template>
        <template #figures>
        </template>
        <!-- FIGURE CAPTION -->
        <template #figureCaption>
        </template>
        <!-- EXPLANATION -->
        <template #belowExplanation>
        </template>
    </VizSection>
</template>

<script setup>
    import { onMounted } from "vue";
    import * as d3 from 'd3';
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        try {
            console.log('app')
        } catch (error) {
            console.error('Error during component mounting', error);
        }        
    });
</script>

<style scoped lang="scss">
</style>
<style lang="scss">
/* css for elements added/classed w/ d3 */
</style>