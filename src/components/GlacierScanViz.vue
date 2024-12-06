<template>
    <section>
        <VizSection
            id="cross-section"
            :figures="true"
            :fig-caption="false"
        >
            <template #figures>
                <div id="cross-section-grid-container">                
                    <div id="caption-container" >
                        <img v-if="defaultView" id="globe-image" src="https://labs.waterdata.usgs.gov/visualizations/images/FireInIce/globe_marker_40.png" alt="locator map showing location of Juneau ice field in Alaska">
                        <img v-if="!defaultView" class="jif-image" :id=currentPhotoID :src=getImageSrc(currentPhotoID) alt="currentPhotoAlt">
                        <div v-if="!mobileView && defaultView">
                            <p v-html="text.paragraph1" />
                            <p v-html="text.promptDesktop" />
                        </div>
                        <div v-if="mobileView && defaultView">
                            <p v-html="text.paragraph1Mobile" />
                            <p v-html="text.promptMobile" />
                        </div>
                        <p v-if="!defaultView" v-html="currentPhotoText"></p>
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
                    {{ text.heading }}
                </h2>
            </template>
            <template #aboveExplanation>
                <p v-html="text.paragraph2" />
                <p v-html="text.paragraph3" />
                <p v-html="text.paragraph4" />
                <p v-html="text.paragraph5" />
            </template>
        </VizSection>
    </section>
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
    const currentXsID = ref("113");
    const currentPhotoID = ref(-9999)
    const currentPhotoAlt = ref("")
    const currentPhotoText = ref("")
    const defaultView = ref(true)
    const default_xs = "113";
    const defaultPhotoID = -9999;

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        try {
            // Use external svg from s3
            d3.xml("https://labs.waterdata.usgs.gov/visualizations/svgs/glacial_mri_v10.svg").then(function(xml) {
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
        currentPhotoAlt.value = props.text[`photo${photo_id}Alt`];
        defaultView.value = false;//!defaultView.value
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
            defaultView.value = true;
        }
    }

    function touchstart(event) {
        draw_xs(default_xs, defaultPhotoID);
        if (event.currentTarget.id.startsWith("figure_1")){
            remove_xs(default_xs, currentPhotoID.value);
            d3.select("#tutorial_arrow").selectAll("path")
                .style("opacity", 0);
        }
        if (event.currentTarget.id.startsWith("xs-main-")){
            remove_xs(currentXsID.value, currentPhotoID.value);
            const line_id = event.currentTarget.id.slice(8);
            currentXsID.value = line_id;
            draw_xs(line_id, defaultPhotoID);
        } else if (event.currentTarget.id.startsWith("xs-c-lg-")){
            remove_xs(currentXsID.value, currentPhotoID.value);
            const line_id = event.currentTarget.id.slice(8);
            currentXsID.value = line_id;
            draw_xs(line_id, defaultPhotoID);
        } else if (event.currentTarget.id.startsWith("photo-lg-")){
            remove_xs(currentXsID.value, currentPhotoID.value);
            const line_id = event.currentTarget.id.slice(13);
            currentXsID.value = line_id;
            const photo_id = event.currentTarget.id.slice(9,12);
            currentPhotoID.value = photo_id;
            draw_xs(line_id, currentPhotoID.value);
            draw_image(photo_id);
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
        if (mobileView) {
            // add event listener to chart groups on mobile to track taps ON svg groups
            // this covers mouseenter, mouseover, and mouseout behavior
            svg.selectAll("g")
                .on("touchstart",(event) => {
                    event.preventDefault();
                    touchstart(event)
                })     
            // add event listener to document to track tap OFF of svg
            // this covers mouseleave behavior
            document.addEventListener('touchstart', function(event) {
                event.preventDefault();
                if (!event.target.ownerSVGElement) {
                    remove_xs(currentXsID.value, currentPhotoID.value);
                    draw_xs(default_xs, defaultPhotoID);
                    d3.select("#tutorial_arrow").selectAll("path")
                        .style("opacity", 0.75);
                    defaultView.value = true;
                }
            }, false);       
        } else {
            svg.selectAll("g")
                .on("mouseover", (event) => mouseover(event))
                .on("mouseout", (event) => mouseout(event))
                .on("mouseenter", (event) => mouseenter(event))
                .on("mouseleave", (event) => mouseleave(event));
        }
    }
</script>

<style scoped lang="scss">
    #cross-section-grid-container {
        display: flex;
        flex-direction: row;
        max-width: 1500px;
        margin: 3rem auto 4rem auto;
        @media screen and (max-width: 600px) {
            flex-direction: column;
        }
    }
    #caption-container {
        height: 100%;
        width: 35vw;
        max-width: 600px;
        margin: auto 3rem auto 0;
        display: flex;
        flex-direction: column;
        flex-grow: 0;
        flex-shrink: 0;
        align-items: center;
        padding: 0rem 3rem 2rem 3rem;
        background-color: var(--faded-usgs-blue);
        border-radius: 5px;
        box-shadow: 5px 5px 10px rgba(57, 61, 66, 0.2);
        font-style: italic;
        @media screen and (max-height: 770px) {
            height: 100%;
            width: 40vw;
            max-width: 40vw;
            padding: 0rem 2rem 1rem 2rem;
        }
        @media screen and (max-width: 600px) {
            flex-direction: column;
            width: 100%;
            max-width: 100%;
            height: 100%;
            padding: 1rem 1.5rem 0 1.5rem;
        }
    }
    .jif-image {
        height: 100%;
        pointer-events: none;
        border-radius: 5px;
        margin-right: 3rem;
        max-width: 100%;
        max-height: 45vh;
        margin: 2rem 2rem 2rem 2rem;
        @media screen and (max-width: 600px) {
            padding-right: 0rem;
            max-height: 35vh;
            max-width: 100%;
            margin: 1rem;
        }
    }
    #globe-image {
        width: 15vw;
        max-width: 250px;
        align-items: center;
        pointer-events: none;
        border-radius: 5px;
        margin: 4rem 2rem 2rem 2rem;
        @media screen and (max-width: 600px) {
            width: 35vw;
            margin: 1rem;
        }
    }
</style>
<style lang="scss">
/* css for elements added/classed w/ d3 */
    #cross-section-svg {
        height: 60vh;
        z-index: 1;
        margin-top: 3rem;
        @media screen and (max-height: 770px) {
            height: 80vh;
        }
        @media screen and (max-width: 600px) {
            height: 45vh;
        }
    } 
</style>