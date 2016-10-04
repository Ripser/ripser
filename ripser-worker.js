importScripts('emscripten/ripser.js');

addEventListener('message', function(e) {
    var data = e.data;
    tic = (new Date()).getTime();
    Module.ripser_emscripten(data.file, data.dim, data.threshold, data.format);
    toc = (new Date()).getTime();

	postMessage(null);
}, false);
