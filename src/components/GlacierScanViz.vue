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
            <div id="cross_section-grid-container">
                <cross_sectionPlot
                    id="cross_section-svg2"
                />
                <svg
                    v-bind="svgAttributes"
                    v-html="svgContent"
                    id="cross_section-svg"
                >
                </svg>
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
    import cross_sectionPlot from "@/assets/svgs/cross_section.svg";
    import VizSection from '@/components/VizSection.vue';

    // define props
    defineProps({
        text: { type: Object }
    })

    // global variables
    const svgName = 'cross_section'
    let svgPromise = ref()
    let svgResult = ref()
    let svgContent = ref()
    let svgAttributes = ref()
    // const {content, attributes} = parseSvg(svgName)

    // console.log(content)
    // console.log(attributes)
    // const importedSvg = ref(null)

    // onBeforeMount(async () => {
    //     importedSvg.value = await parseSvg('cross_section')
    // })
    // console.log(importedSvg)

    // Declare behavior on mounted
    // functions called here
    onMounted(async () => {
        try {
            await parseSvg(svgName)
            
            svgPromise.value.then(function(result) {
                addSvgContent(result);
                addInteractions();
            })
        } catch (error) {
            console.error('Error during component mounting', error);
        }        
    });

    // load svg from s3
    async function parseSvg(svgString) {
        let url = `https://labs.waterdata.usgs.gov/visualizations/svgs/${svgString}.svg`
        try {
            svgPromise.value = fetch(url)
                .then((response) => response.text())
                .then((svgText) => {
                    const parsedSvg = new DOMParser()
                        .parseFromString(svgText, 'image/svg+xml')
                        .querySelector('svg')
                    
                    return(parsedSvg)
                })
            console.log('svg parsed')
        } catch (error) {
            console.error('Error parsing svg', error);
        }

        

        // const svgElement = fetch(url)
        //     .then((response) => response.text())
        //     .then((svgText) => {
        //         const parsedSvg = new DOMParser()
        //             .parseFromString(svgText, 'image/svg+xml')
        //             .querySelector('svg')
                
        //         return(parsedSvg)
        //     })

        // return(svgElement)


        // return {
        //     content: svgElement.innerHTML,
        //     attributes: Object.fromEntries([...svgElement.attributes].map(a => [a.name, a.value]))
        // }
    }

    function addSvgContent(result) {      
        svgResult.value = result;
        // console.log(svgResult.value.innerHTML)
        // console.log('that was add svg content')

        svgContent.value = svgResult.value.innerHTML.replaceAll(' xmlns="http://www.w3.org/2000/svg"','');
        // console.log(svgContent.value)
        // console.log(svgContent.value.replaceAll(' xmlns="http://www.w3.org/2000/svg"',''))
        svgAttributes.value = Object.fromEntries([...svgResult.value.attributes].map(a => [a.name, a.value]))
        console.log('Added svg content')
        // return {
        //     content: svgResult.value.innerHTML,
        //     attributes: Object.fromEntries([...svgResult.value.attributes].map(a => [a.name, a.value]))
        // }   
    }

    function draw_xs(line_id){
        console.log(line_id)
        console.log(d3.select("#xs-main-" + line_id).selectAll("path"))
        console.log(d3.select("#xs-topo-" + line_id).selectAll("path"))
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
            .style("fill", "#D2B48C")
            .style("fill-opacity", 1);
        d3.select("#xs-ice-" + line_id).selectAll("path")
            .style("fill", "#b9e8ea")
            .style("fill-opacity", 1);
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
            .style("fill-opacity", 0);
        d3.select("#xs-ice-" + line_id).selectAll("path")
            .style("fill-opacity", 0);
    }   

    function mouseover(event,default_xs) {
        if (event.currentTarget.id.startsWith("xs-main-")){
            remove_xs(default_xs)
            let line_id = event.currentTarget.id.slice(8);
            draw_xs(line_id)
        }
    }

    function mouseout(event) {
        if (event.currentTarget.id.startsWith("xs-main-")){
            let line_id = event.currentTarget.id.slice(8);
            remove_xs(line_id)
        }
    }

    function mouseleave(event,default_xs) {
        if (event.currentTarget.id.startsWith("figure_1")){
            draw_xs(default_xs)
        }
    }

    function addInteractions() {
        // set viewbox for svg with loss function chart
        const cross_sectionSVG = d3.select("#cross_section-svg");
        console.log('whole svg:')
        console.log(document.getElementById("cross_section-svg2"))
        console.log(cross_sectionSVG)
        console.log('orig svg:')
        console.log(d3.select("#cross_section-svg2"))
        console.log('svg groups:')
        console.log(cross_sectionSVG.selectAll("g"))
        console.log('svg selection')
        console.log(cross_sectionSVG.select("#xs-main-50").selectAll("path"))

        var default_xs = "50"
        draw_xs(default_xs)

        // // Add interaction to loss function chart
        // cross_sectionSVG.selectAll("g")
        //     .on("mouseover", (event) => mouseover(event,default_xs))
        //     .on("mouseout", (event) => mouseout(event))
        //     .on("mouseleave", (event) => mouseleave(event,default_xs))
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
    #cross_section-svg {
        grid-area: chart;
        place-self: center;
        max-height: 80%;
        max-width: 80%;
    }

</style>