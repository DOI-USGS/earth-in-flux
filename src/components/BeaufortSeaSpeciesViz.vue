<template>
    <!---VizSection-->
    <VizSection
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
        </template>
        <template #figures>
            <div class="single maxWidth">
                <figure>
                    <img src="../../beaufortSea/out/BeaufortSeaSpecies_1.png">
                </figure>
                <p>As climate conditions change at the mouth of the Mackenzie River, where it spills into the Beaufort Sea north of Alaska, the small species of ostracodes and forams respond and the overall assemblage of species changes. Researchers have reconstructed the past climates of the Arctic using these microfossils. In particular, a couple species work as indicator species, which means we can monitor their growth (or lack thereof) to understand how healthy the local ocean climate is at a given time in the past or present. </p>
                <h1>Meet the ostracodes!</h1>
                <figure>
                    <img src="../../beaufortSea/out/BeaufortSeaSpecies_2.png">
                </figure>
                <p>Ostracodes are a microscopic species of arthropods, distantly related to crustaceans and insects. The fossil remains are from their hardened shells, each distinctly patterned. Kotorachythere arctoborealis is a species that was relatively sparse in the Beaufort Sea fossil record until more recent years. Paracyprideis pseudopunctillata is one of the most dominant species in this Arctic fossil record, at least until the recent decades.</p>
                <h1>Meet the forams!</h1>
                <figure>
                    <img src="../../beaufortSea/out/BeaufortSeaSpecies_3.png">
                </figure>
                <p>Foraminifera, or "forams", are a group of single-celled protists covered by an external shell called a "test." It is these tests that are found in the sedimentary fossil record and are used to identify the different species. The tests, like those of the Cassidulina and Elphidium species, are often made of calcium carbonate. However, some species such as Spiroplectammina create tests that are agglutinated, or clumped together particles such as sand and sediment.</p>
            </div>
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
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })
</script>

<style>
</style>