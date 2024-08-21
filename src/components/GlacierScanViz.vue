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
            <div id="cross_section-grid-container" />
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
            // Use external svg from s3
            d3.xml("https://labs.waterdata.usgs.gov/visualizations/svgs/glacial_xray.svg").then(function(xml) {
                // add svg content to DOM
                const svgGrid = document.getElementById("cross_section-grid-container")
                svgGrid.appendChild(xml.documentElement);
                
                // add id to svg
                d3.select("#cross_section-grid-container").select("svg")
                    .attr("id", "cross_section-svg")

                // add interactivity
                addInteractions();
            });
        } catch (error) {
            console.error('Error during component mounting', error);
        }        
    });

    function draw_xs(line_id){
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
            .style("fill", "#4fd437")
            .style("fill-opacity", 1)
            .style("stroke", "#000000")
            .style("stroke-opacity", 1.0);
    }
    
    function remove_xs(line_id){
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
    }   

    function mouseover(event,default_xs) {
        if (event.currentTarget.id.startsWith("xs-main-")){
            remove_xs(default_xs);
            let line_id = event.currentTarget.id.slice(8);
            draw_xs(line_id);
        } else if (event.currentTarget.id.startsWith("xs-c-lg-")){
            remove_xs(default_xs);
            let line_id = event.currentTarget.id.slice(8);
            draw_xs(line_id);
        }
    }

    function mouseout(event) {
        if (event.currentTarget.id.startsWith("xs-main-")){
            let line_id = event.currentTarget.id.slice(8);
            remove_xs(line_id);
        } else if (event.currentTarget.id.startsWith("xs-c-lg-")){
            let line_id = event.currentTarget.id.slice(8);
            remove_xs(line_id);
        }
    }

    function mouseenter(event,default_xs) {
        if (event.currentTarget.id.startsWith("figure_1")){
            remove_xs(default_xs);
            d3.select("#tutorial_arrow").selectAll("path")
                .style("opacity", 0);
        }
    }

    function mouseleave(event,default_xs) {
        if (event.currentTarget.id.startsWith("figure_1")){
            draw_xs(default_xs);
            d3.select("#tutorial_arrow").selectAll("path")
                .style("opacity", 0.75);
        }
    }

    function addInteractions() {
        // set viewbox for svg with loss function chart
        const cross_sectionSVG = d3.select("#cross_section-svg");

        var default_xs = "188";
        draw_xs(default_xs);

        // Add interaction to loss function chart
        cross_sectionSVG.selectAll("g")
            .on("mouseover", (event) => mouseover(event,default_xs))
            .on("mouseout", (event) => mouseout(event))
            .on("mouseenter", (event) => mouseenter(event,default_xs))
            .on("mouseleave", (event) => mouseleave(event,default_xs));
    }
</script>

<style scoped lang="scss">
    #cross_section-grid-container {
        display: grid;
        width: 100%;
        max-width: 1200px;
        margin: 0 auto 0 auto;
        grid-template-areas:
            "chart";
    }
</style>
<style lang="scss">
/* css for elements added/classed w/ d3 */
    #cross_section-svg {
        grid-area: chart;
        place-self: center;
        max-height: 100%;
        max-width: 100%;
    }
</style>