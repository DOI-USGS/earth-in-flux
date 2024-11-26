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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_1.png" alt="Illustrated world map showing a cartoon airplane flying from Reston, Virginia at USGS headquarters to the Beaufort Sea at the northeast corner of Alaska.">
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_2.png" alt="Zoomed in illustrated map of Alaska and northwest Canada. A cartoon icebreaker ship is offshore from the mouth of the Mackenzie River where it opens up into the Beaufort Sea. Text labels Alaska, Canada, Mackenzie River, and Beaufort Sea.">
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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_3.gif" alt="Animated cartoon gif showing the icebreaker ship pushing through some ice to get to the study site. The landscape shows ice-covered Beaufort Sea with a dull, cold blue-grey sky and distant clouds.">
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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_4.png" alt="Two side-by-side photographs of the ship on the left and people on deck of the ship on the right. The images are composed in a scrapbook-style alongside a cartoon beluga whale. The ship is large and bright red with the U.S. Coast Guard label on the side. Behind the ship are mountains covered in ice. The people on the deck of the ship are watching the sediment core as it enters the ice. They're wearing bright red Coast Guard jackets to stay safe and warm.">
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_5.gif" alt="Animated cartoon gif showing the Coast Guard ship on top of the icy sea and dropping the sediment core down into the ocean. When the core gets to the bottom, it removes a core of the ocean floor and then is returned to the ship with the sediment sample. Swimming under the sea is a cartoon beluga whale.">
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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_6.png" alt="Three photographs of the researchers collecting the core, processing it, and then also showing the core once it is in the lab. The images are composed in a scrapbook-style alongside the cartoon Coast Guard ship.">
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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_7.gif" alt="Animated cartoon gif showing the sediment core after it's been pulled up from the ocean floor. As the animation moves forward, the core is sliced into small discs. One example disc is shown zoomed-in with lots of little specks of color, representing the microfossils. The five focal species are also shown zoomed-in near the sample to represent that they are identified from these sediment discs.">
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
                    <img class="coring-image" src="https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_8.png" alt="Three photographs of the researchers examining and identifying the microfossils. One image shows Laura Gemery, USGS researcher, studying microfossils through a microscope. The middle image shows a zoomed in image of a petri dish and Laura holding a small paint brush and one microfossil, which is barely visible next to the paint brush bristles. The third image shows a close-up of the microscope with some samples laid out underneath it. The scrapbook-style images are surrounded by actual scanning electron images of the different microfossil species.">
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
    max-width: 70rem;
    justify-self: center;
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