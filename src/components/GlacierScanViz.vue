<template>
    <VizSection
        id="cross-section"
        :figures="true"
        :fig-caption="false"
    >
        <template #heading>
            <h2>
                {{ text.heading1 }}
            </h2>
        </template>
        <template #figures>
            <div id="cross-section-grid-container">                
                <div id="caption-container">
                    <p v-html="text.paragraph1" />
                    <p v-if="!mobileView" v-html="text.promptDesktop" />
                    <p v-if="mobileView" v-html="text.promptMobile" />
                </div>
                <div id="note-container" class="hide">
                    <p v-html="currentPhotoText"></p>
                </div>
                <div id="photo-container" class="hide">
                    <img class="jif-image" :id=currentPhotoID :src=getImageSrc(currentPhotoID) alt="currentPhotoAlt">
                </div>
            </div>
        </template>
    </VizSection>

    <VizSection
        id="cross-section-how-to"
        :figures="false"
        :fig-caption="false"
    >
        <template #heading>
            <h2>
                {{ text.heading2 }}
            </h2>
        </template>
        <template #aboveExplanation>
            <p v-html="text.paragraph2" />
            <p v-html="text.paragraph3" />
        </template>
    </VizSection>
</template>

<script setup>
    import { onMounted, ref } from "vue";
    import * as d3 from 'd3';
    import { isMobile } from 'mobile-device-detect';
    import VizSection from '@/components/VizSection.vue';

    // define props
    const props = defineProps({
        text: { type: Object }
    })

    // global variables
    const mobileView = isMobile;
    const currentPhotoID = ref(null)
    const currentPhotoAlt = ref("")
    const currentPhotoText = ref("")
    const default_xs = "113";

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        try {
            // Use external svg from s3
            d3.xml("https://labs.waterdata.usgs.gov/visualizations/svgs/glacial_mri_v9.svg").then(function(xml) {
                // add svg content to DOM
                const svgGrid = document.getElementById("cross-section-grid-container")
                svgGrid.appendChild(xml.documentElement);

                // add id to svg
                const crossSectionSVG = d3.select("#cross-section-grid-container").select("svg")
                    .attr("id", "cross-section-svg")
                    .attr("width", "100%")
                    .attr("height", "100%")

                // hide some components
                crossSectionSVG.select("#tutorial-dt-1")
                    .attr("display", "none")
                crossSectionSVG.select("#tutorial-dt-2")
                    .attr("display", "none")
                crossSectionSVG.select("#tutorial-mb-1")
                    .attr("display", "none")
                crossSectionSVG.select("#tutorial-mb-2")
                    .attr("display", "none")
                crossSectionSVG.select("#legend_1")
                    .style("transform", "translate(-100px, 210px)")
                crossSectionSVG.select("#legend_1").select("#patch_3")
                    .attr("display", "none")

                // add interactivity
                addInteractions(crossSectionSVG);
            });
        } catch (error) {
            console.error('Error during component mounting', error);
        }        
    });

    function getImageSrc(photoID) {
        return `https://labs.waterdata.usgs.gov/visualizations/images/FireInIce/juneau_icefield_${photoID}.jpeg`
    }

    function draw_xs(line_id,photo_id){
        d3.select("#xs-main-" + line_id).selectAll("path")
            .style("stroke-opacity", 1.0);
        d3.select("#xs-w" + line_id).selectAll("path")
            .style("stroke-opacity", 1.0);
        d3.select("#xs-label1-"+ line_id).selectAll("text")
            .style("opacity", 1)
            .style("font-weight", 10);
        d3.select("#xs-w-label1-"+ line_id).selectAll("text")
            .style("opacity", 1)
            .style("font-weight", 1000);
        d3.select("#xs-arrow1-"+ line_id).selectAll("path")
            .style("fill-opacity", 1);
        d3.select("#xs-w-arrow1-"+ line_id).selectAll("path")
            .style("fill-opacity", 1);
        d3.select("#xs-label2-"+ line_id).selectAll("text")
            .style("opacity", 1);
        d3.select("#xs-w-label2-"+ line_id).selectAll("text")
            .style("opacity", 1);
        d3.select("#xs-arrow2-"+ line_id).selectAll("path")
            .style("fill-opacity", 1);
        d3.select("#xs-w-arrow2-"+ line_id).selectAll("path")
            .style("fill-opacity", 1);
        d3.select("#xs-topo-" + line_id).selectAll("path")
            .style("fill", "#c49051")
            .style("fill-opacity", 1)
            .style("stroke", "#000000")
            .style("stroke-opacity", 0.1);
        d3.select("#xs-ice-" + line_id).selectAll("path")
            .style("fill", "#dddddd")
            .style("fill-opacity", 1)
            .style("stroke", "#000000")
            .style("stroke-opacity", 0.2);
        d3.select("#xs-c-sm-" + line_id).selectAll("path")
            .style("fill", "#1f77b4")
            .style("fill-opacity", 1)
            .style("stroke", "#000000")
            .style("stroke-opacity", 1.0);
        d3.select("#photo-sm-"+photo_id +"-"+ line_id).selectAll("path")
            .style("fill", "#d62728")
            .style("fill-opacity", 0.8)
            .style("stroke", "#000000")
            .style("stroke-opacity", 1.0);
    }
    
    function remove_xs(line_id,photo_id){
        d3.select("#xs-main-" + line_id).selectAll("path")
            .style("stroke-opacity", 0);
        d3.select("#xs-w" + line_id).selectAll("path")
            .style("stroke-opacity", 0);
        d3.select("#xs-label1-"+ line_id).selectAll("text")
            .style("opacity", 0);
        d3.select("#xs-w-label1-"+ line_id).selectAll("text")
            .style("opacity", 0);
        d3.select("#xs-arrow1-"+ line_id).selectAll("path")
            .style("fill-opacity", 0);
        d3.select("#xs-w-arrow1-"+ line_id).selectAll("path")
            .style("fill-opacity", 0);
        d3.select("#xs-label2-"+ line_id).selectAll("text")
            .style("opacity", 0);
        d3.select("#xs-w-label2-"+ line_id).selectAll("text")
            .style("opacity", 0);
        d3.select("#xs-w-arrow2-"+ line_id).selectAll("path")
            .style("fill-opacity", 0);
        d3.select("#xs-arrow2-"+ line_id).selectAll("path")
            .style("fill-opacity", 0);
        d3.select("#xs-topo-" + line_id).selectAll("path")
            .style("fill-opacity", 0)
            .style("stroke-opacity", 0);
        d3.select("#xs-ice-" + line_id).selectAll("path")
            .style("fill-opacity", 0)
            .style("stroke-opacity", 0);
        d3.select("#xs-c-sm-" + line_id).selectAll("path")
            .style("fill-opacity", 0)
            .style("stroke-opacity", 0);
        d3.select("#photo-sm-"+photo_id +"-"+ line_id).selectAll("path")
            .style("fill-opacity", 0)
            .style("stroke-opacity", 0);
    }   

    function draw_image(photo_id){
        currentPhotoID.value = photo_id
        currentPhotoText.value = mobileView ? props.text[`photo${photo_id}Mobile`] : props.text[`photo${photo_id}`];
        d3.select("#caption-container").selectAll("p")
            .attr("class", "hide");
        d3.select("#photo-container")
            .attr("class", "show");
        d3.select("#note-container")
            .attr("class", "show");
    }

    function remove_image(){
        d3.select("#caption-container").selectAll("p")
            .attr("class", "show");
        d3.select("#photo-container")
            .attr("class", "hide");
        d3.select("#note-container")
            .attr("class", "hide");
    }

    function mouseover(event) {
        if (event.currentTarget.id.startsWith("xs-main-")){
            remove_xs(default_xs,-9999);
            const line_id = event.currentTarget.id.slice(8);
            draw_xs(line_id,-9999);
        } else if (event.currentTarget.id.startsWith("xs-c-lg-")){
            remove_xs(default_xs,-9999);
            const line_id = event.currentTarget.id.slice(8);
            draw_xs(line_id,-9999);
        } else if (event.currentTarget.id.startsWith("photo-lg-")){
            remove_xs(default_xs,-9999);
            const line_id = event.currentTarget.id.slice(13);
            const photo_id = event.currentTarget.id.slice(9,12);
            draw_xs(line_id,photo_id);
            draw_image(photo_id);
        }
    }

    function mouseout(event) {
        if (event.currentTarget.id.startsWith("xs-main-")){
            const line_id = event.currentTarget.id.slice(8);
            remove_xs(line_id,-9999);
        } else if (event.currentTarget.id.startsWith("xs-c-lg-")){
            const line_id = event.currentTarget.id.slice(8);
            remove_xs(line_id,-9999);
        } else if (event.currentTarget.id.startsWith("photo-lg-")){
            const line_id = event.currentTarget.id.slice(13);
            const photo_id = event.currentTarget.id.slice(9,12);
            remove_xs(line_id,photo_id);
            remove_image();
        }
    }

    function mouseenter(event) {
        if (event.currentTarget.id.startsWith("figure_1")){
            remove_xs(default_xs,-9999);
            d3.select("#tutorial_arrow").selectAll("path")
                .style("opacity", 0);
        }
    }

    function mouseleave(event) {
        if (event.currentTarget.id.startsWith("figure_1")){
            draw_xs(default_xs,-9999);
            d3.select("#tutorial_arrow").selectAll("path")
                .style("opacity", 0.75);
        }
    }

    function addInteractions(svg) {

        if (mobileView == true){
            d3.select('#tutorial-dt-1').selectAll("path")
                .style("opacity", 0.0);
            d3.select('#tutorial-dt-1').selectAll("text")
                .style("opacity", 0.0);
            d3.select('#tutorial-dt-2').selectAll("path")
                .style("opacity", 0.0);
            d3.select('#tutorial-dt-2').selectAll("text")
                .style("opacity", 0.0);
            d3.select('#tutorial-mb-1').selectAll("path")
                .style("opacity", 0.75);
            d3.select('#tutorial-mb-1').selectAll("text")
                .style("opacity", 1.0);
            d3.select('#tutorial-mb-2').selectAll("path")
                .style("opacity", 0.75);
            d3.select('#tutorial-mb-2').selectAll("text")
                .style("opacity", 1.0);
        } else{
            d3.select('#tutorial-dt-1').selectAll("path")
                .style("opacity", 0.75);
            d3.select('#tutorial-dt-1').selectAll("text")
                .style("opacity", 1.0);
            d3.select('#tutorial-dt-2').selectAll("path")
                .style("opacity", 0.75);
            d3.select('#tutorial-dt-2').selectAll("text")
                .style("opacity", 1.0);
            d3.select('#tutorial-mb-1').selectAll("path")
                .style("opacity", 0.0);
            d3.select('#tutorial-mb-1').selectAll("text")
                .style("opacity", 0.0);
            d3.select('#tutorial-mb-2').selectAll("path")
                .style("opacity", 0.0);
            d3.select('#tutorial-mb-2').selectAll("text")
                .style("opacity", 0.0);
        }

        draw_xs(default_xs,-9999);

        // Add interaction events
        svg.selectAll("g")
            .on("mouseover", (event) => mouseover(event))
            .on("mouseout", (event) => mouseout(event))
            .on("mouseenter", (event) => mouseenter(event))
            .on("mouseleave", (event) => mouseleave(event));
    }
</script>

<style scoped lang="scss">
    #cross-section-grid-container {
        display: grid;
        width: 100%;
        max-width: 1200px;
        max-height: 85vh;
        margin: 4rem auto 4rem auto;
        grid-template-columns: max-content 1fr;
        grid-template-rows: 20vh 1fr;
        row-gap: 2rem;
        grid-template-areas:
            "photo text"
            "chart chart";
        @media screen and (max-height: 770px) {
            max-height: 80vh;
            column-gap: 5rem;
            grid-template-columns: 50% 50%;
            grid-template-rows: 60% 40%;
            grid-template-areas:
                "photo chart"
                "text chart";
        }
        @media screen and (max-width: 600px) {
            max-height: 100vh;
            grid-template-columns: 100%;
            grid-template-rows: 20% 15% 65%;
            row-gap: 1rem;
            grid-template-areas:
                "photo"
                "text"
                "chart";
        }
    }
    #caption-container {
        grid-column: 1 / span 2;
        grid-row: 1;
        background-color: var(--faded-usgs-blue);
        border-radius: 5px;
        box-shadow: 5px 5px 10px rgba(57, 61, 66, 0.2);
        align-content: center;
        font-style: italic;
        padding: 0 5rem 0 5rem;
        @media screen and (max-height: 770px) {
            grid-column: 1;
            grid-row: 1 / span 2;
        }
        @media screen and (max-width: 600px) {
            grid-column: 1;
            grid-row: 1 / span 2;
            padding: 1rem 1rem 0 1rem;
        }
    }
    #note-container {
        grid-area: text;
        align-content: center;
        font-style: italic;
        padding: 2rem;
        @media screen and (max-height: 770px) {
            padding: 0 4rem 2rem 4rem;
            align-content: start;
        }
        @media screen and (max-width: 600px) {
            padding: 0 1rem 1rem 1rem;
            align-content: start;
        }
    }
    #note-container p {
        padding: 0rem;
    }
    #photo-container {
        grid-area: photo;
        align-content: center;
        text-align: center;
        border-radius: 5px;
        @media screen and (max-height: 770px) {
            padding-top: 2rem;
        }
        @media screen and (max-width: 600px) {
            padding-top: 1rem;
        }
    }
    .jif-image {
        align-items: center;
        max-height: 100%;
        pointer-events: none;
        border-radius: 5px;
    }
    .hide {
        opacity: 0;
        transition: opacity ease-out 0.2s;
    }
    .show {
        opacity: 1;
        transition: opacity ease-in 0.2s;
    }
</style>
<style lang="scss">
/* css for elements added/classed w/ d3 */
    #cross-section-svg {
        grid-area: chart;
        place-self: center;
        z-index: 1;
    } 
</style>