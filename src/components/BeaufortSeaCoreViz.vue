<template>
    <div>
        <VizSection
            :figures="true"
            :fig-caption="false"
        >
            <template #heading>
                <h2 v-html="text.heading1" />
            </template>
            <template #aboveExplanation>
                <p v-html="text.intro" />
            </template>
            <template #belowExplanation>
                <h2 v-html="text.heading2" />
            </template>
        </VizSection>
        <div id="sediment-coring-grid-container">
            <div id="coring-image-container">
                <img class="coring-image" :src="currentImage" :alt="currentAltText">
            </div>
            <div id="coring-text-container" class="text-container">
                <p v-html="currentText" />
            </div>
            <button id="coring-prev" class="flip-button" @click="currentIndex--" :disabled="isFirstImage">
                <font-awesome-icon :icon="{ prefix: 'fas', iconName: 'arrow-left' }"  class="fa fa-arrow-left"/>
            </button>
            <button id="coring-next" class="flip-button" @click="currentIndex++" :disabled="isLastImage">
                <font-awesome-icon :icon="{ prefix: 'fas', iconName: 'arrow-right' }"  class="fa fa-arrow-right"/>
            </button>
        </div>
    </div>
</template>

<script setup>
    import { computed, ref } from 'vue';
    import VizSection from '@/components/VizSection.vue';

    // define props
    const props = defineProps({
        text: { type: Object }
    })

    // global variables
    const currentIndex = ref(1);
    const gifImageIndices = [3, 5, 7]
    const nImages = 8;

    const isFirstImage = computed(() => {
        return currentIndex.value === 1;
    });
    const isLastImage = computed(() => {
        return currentIndex.value === nImages;
    });
    const currentText = computed(() => {
        const selectionString = 'paragraph' + currentIndex.value
        return props.text[selectionString];
    });
    const currentImage = computed(() => {
        let fileEnding;
        if (gifImageIndices.includes(currentIndex.value)) {
            fileEnding = 'gif';
        } else {
            fileEnding = 'png'
        }
        return `https://labs.waterdata.usgs.gov/visualizations/images/BeaufortSea/BeaufortSeaCore_${currentIndex.value}.${fileEnding}`;
    });
    const currentAltText = computed(() => {
        const selectionString = 'alt' + currentIndex.value
        return props.text[selectionString];
    });
</script>

<style>
    #sediment-coring-grid-container{
        display: grid;
        max-width: 1200px;
        grid-template-columns: 10% calc(80% - 4rem) 10%;
        grid-template-rows: auto max-content;
        grid-template-areas:
            "prev image next"
            "text text text";
        margin: 2rem auto 0 auto;
        column-gap: 2rem;
        row-gap: 3rem;
        @media only screen and (max-width: 600px) {
            width: 90vw;
            grid-template-rows: auto max-content;
            grid-template-areas:
                "image image image"
                "prev text next";
        }
    }
    #coring-image-container {
        grid-area: image;
        justify-self: center;
        align-self: center;
    }
    .coring-image {
        width: 100%;
        justify-self: center;
        @media only screen and (max-width: 600px) {
            width: 100%;
        }
    }
    #coring-text-container {
        grid-area: text;
        height: 15vh;
        @media screen and (max-height: 770px) {
            height: 30vh;
        }
        @media only screen and (max-width: 600px) {
            height: auto;
        }
    }
    .flip-button {
        height: 5rem;
        width: 5rem;
        align-self: center;
        border-radius: 5rem;
        border: solid 0.75px var(--medium-grey);
        cursor: pointer;
        box-shadow: 0px 0px 4px rgba(39,44,49,.3);
        @media only screen and (max-width: 600px) {
            height: 3rem;
            width: 3rem;
            align-self: start;
        }
    }
    #coring-prev {
        grid-area: prev;
        justify-self: end;
    }
    #coring-next {
        grid-area: next;
        justify-self: start;
    }
    button:hover:after {
        top: 0px;
        left: 0px;
    }
    button:hover {
        box-shadow: rgba(39,44,49,.7) 2px 2px 4px -2px;
        transform: translate3d(0, 2px, 0);
    }
    button:disabled {
        background-color: #b4b2b2;
        cursor: not-allowed;
        border-color: #b4b2b2;
        color: #b4b2b2;
    }
    button:disabled:hover {
        box-shadow: 0px 0px 4px rgba(39,44,49,.3);
        transform: translate3d(0, 0px, 0);
    }
    button:disabled:after {
        background-color: #b4b2b2;
    }
    .fa {
        color: var(--color-text);
        opacity: 1;
        @media only screen and (max-width: 600px) {
            font-size: 1.6rem;
        }
    }
    button:disabled .fa {
        opacity: 0.2;
    }
</style>