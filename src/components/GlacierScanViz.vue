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
        </template>
        <template #figures>
            <div id="cross-section-grid-container">
                <div id="photo-container">
                    <img class="jif-image" :id=currentPhotoID :src=getImageSrc(currentPhotoID) alt="currentPhotoAlt">
                </div>
                <div id="caption-container">
                    <p v-html="currentPhotoText"></p>
                </div>
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
                d3.select("#cross-section-grid-container").select("svg")
                    .attr("id", "cross_section-svg")

                // add interactivity
                addInteractions();
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
        d3.select("#image-container")
            .style("visibility", "visible");
        currentPhotoID.value = photo_id
        currentPhotoText.value = props.text[`photo${photo_id}`]
    }

    function remove_image(){
        d3.select("#image-container")
            .style("visibility", "hidden");
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

    function addInteractions() {
        // set viewbox for svg with loss function chart
        const cross_sectionSVG = d3.select("#cross_section-svg")
            .attr("width", "100%")
            .attr("height", "100%");

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

        // Add interaction to loss function chart
        cross_sectionSVG.selectAll("g")
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
        margin: 4rem auto 0 auto;
        grid-template-columns: 20% 1fr;
        grid-template-rows: 20% 80%;
        grid-template-areas:
            "photo text"
            "chart chart";
       row-gap: 3rem;
       column-gap: 3rem;
        @media screen and (max-height: 770px) {
            grid-template-columns: 30% 70%;
            grid-template-areas:
                "photo chart"
                "text chart";
        }
        @media screen and (max-width: 600px) {
            grid-template-rows: 1fr 1fr 1fr;
            grid-template-columns: 100%;
            grid-template-areas:
                "photo"
                "text"
                "chart";
        } 
    }
    
    #photo-container {
        grid-area: photo;
        
        justify-content: center;
        align-content: start;
    }
    .jif-image {
        max-height: 100%;
        max-width: 100%;
        pointer-events: none;
        border: 2px solid black;
        border-radius: 15px;
        box-shadow: 5px 5px 15px rgba(0, 0, 0, 0.5);
        z-index: 2;
    }
    #caption-container {
        grid-area: text;
    }
</style>
<style lang="scss">
/* css for elements added/classed w/ d3 */
    #cross_section-svg {
        grid-area: chart;
        place-self: center;
        z-index: 1;
    } 
</style>