<template>
    <div class="chart">
        <div class="takeaway">
            <p>{{ description }}</p>
        </div>
        <img :src="src" :alt="alt"/>        
        <div class="project">
            <p>{{ project }}</p>
        </div>

        <!--div class="overlay">
            <p>{{ description }}</p>
        </div-->
    </div>
</template>

<script setup>
    defineProps({
        src: {
            type: String,
            default: ``
        },
        alt: {
            type: String,
            default: ``
        },
        project: {
            type: String,
            default: ''
        },
        description: {
            type: String,
            default: ''
        }
    })

</script>

<style scoped lang="scss">
	.chart{
		display: flex;
		flex-direction: column-reverse;
		justify-content: flex-start;
		position: relative;
		width: 280px;
		height: 28rem;//100%;
        -moz-transition: height 1s ease;
        -webkit-transition: height 1s ease;
        -o-transition: height 1s ease;
        transition: height 1s ease;
        // overflow: hidden;
        background-color: white;
        // align-items: start;
	}
    .chart:hover {
        // transform: scale(1.1);
        height: 370px;
        transform: scaleY(1);
        // overflow: auto;
    }

    .overlay {
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background-color: rgba(0, 0, 0, 0.5);
        color: white;
        opacity: 0;
        display: flex;
        font-weight: 700;
    }

    .chart:hover .overlay {
        opacity: 1;
    }
    .takeaway {
        position: absolute;
        top: 230px;
        left: 0;
        right: 0;
        bottom: 0;
        line-height: 1;
        display: inline-block;
        vertical-align: middle;
        -webkit-transition: top 1s ease-out; 
        -moz-transition: top 1s ease-out; 
        -ms-transition: top 1s ease-out; 
        -o-transition: top 1s ease-out; 
        transition: top 1s ease-out;
    }
    .takeaway p {
        padding: 0.2rem 2rem 0.2rem 2rem; 
    }
    .chart:hover .takeaway {
        top: 300px;
    }
    .chart img {
        position: absolute;
        width: 100%;
        top: 0;
    }

    .overlay p {
        margin: 0;
        padding: 10px;
        text-align: center;
    }
    .project {
        position: absolute;
        top: calc(280px - 2rem);
        width: auto;
        // left: 50%;
    }
    .project p {
        background-color: var(--medium-grey);
        padding: 0.2rem 0.5rem 0.2rem 0.5rem; 
        text-align: center;
        width: auto;
        font-weight: 200;
        border: solid var(--color-background);
        border-radius: 10% 10% 0% 0%;
        // color: var(--color-background)
    }
</style>