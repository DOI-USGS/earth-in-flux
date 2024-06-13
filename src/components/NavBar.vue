<template>
    <header id="NavBar" class="header-nav"  role="banner">
        <div class="topnav">
            <a class="active" href="/">Home</a>
            <a href="/">Takeaways</a>
            <a href="/">Methods</a>
            <a href="/">About the Projects</a>
            <a href="https://labs.waterdata.usgs.gov/visualizations/index.html#/" target="_blank">About VizLab</a>
        </div>
    </header>
</template>

<script setup>
</script>

<style scoped>

.header-nav {
    position: -webkit-sticky; /* For Safari */
    position: sticky;
    top: 0;
    z-index: 1000;
}

.topnav {
    background-color: #6b7c94;
    overflow: hidden;
    display: flex;
    justify-content: center;
}

/* Style the links inside the navigation bar */
.topnav a {
    color: #f2f2f2;
    text-align: center;
    padding: 14px 16px;
    text-decoration: none;
    font-size: 17px;

}

/* Change the color of links on hover */
.topnav a:hover {
    background-color: #d2e9ff;
    color: var(--dodger-blue);
}

/* Add a color to the active/current link */
.topnav a.active {
    background-color: var(--color-title-text);
    color: white;
}
</style>
