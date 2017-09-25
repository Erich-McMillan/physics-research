function providePlugin(pluginName, pluginConstructor) {
    var ga = window[window['GoogleAnalyticsObject'] || 'ga'];
    if (ga) ga('provide', pluginName, pluginConstructor);
}

var googleAnalyticsPlugin = function (tracker, config) {

};

googleAnalyticsPlugin.prototype.doStuff = function () {

};


//googleAnalyticsPlugin.prototype.doStuff2 = function (enabled) {
//alert('doStuff2');
//};

providePlugin('googleAnalyticsPlugin', googleAnalyticsPlugin);