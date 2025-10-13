
function ToggleErrorDisplay(selector) {
    //$('#refTab').toggleClass('hide-error');
    const target = document.querySelector(selector);
    _toggleClasse(target,'hide-error');
}

const _base =[];
function _each (target, cb){
    return _base.forEach.call(target, cb);
}
function _toggleClasse(target, className){
    const rClass = new RegExp (' '+ className + ' ');
    let targetClass = ' ' + target.className +' ';
    if(rClass.test(targetClass)){
        targetClass = targetClass.replace(rClass, ' ' );
    }
    else{
        targetClass += className;
    }
    return target.className = targetClass.trim();
};

function _addClasse(target, className){
    const rClass = new RegExp (' '+ className + ' ');
    let targetClass = ' ' + target.className +' ';
    if(!rClass.test(targetClass)){
        targetClass += className;
    }
    
    return target.className = targetClass.trim();
};

function makeTree(selector){

    const collection = {
        root:{                        
            children: []
        }
    };
    const table = document.querySelector(selector);
    const rTestError = /_error/;
    const split = /(.+)\.([^\.]+)$/;
    const slash = /\//g;
    let tBody;

    _addClasse(table, 'hide-error');

    for(let i = 0; i<table.children.length; i++){
        tBody =table.children[i];
        if (tBody.tagName === "TBODY") break;
    }

    _each(tBody.children , function(row){
    
        
        const pathElem = row.firstElementChild.getElementsByClassName('pathname');
        const fullPath = pathElem[0].textContent.replace(slash, '.');
        const res = fullPath.match(split);
        const entry = collection[fullPath]= {
        
        };
        row.setAttribute("data-tt-id",fullPath);
        // si res est null, il n'y a pas de parent
        if (res === null){
            collection.root.children.push(entry);
        }
        else{
        	pathElem[0].textContent = res[2];
            //row.firstElementChild.textContent = res[2];
            const parent = collection[res[1]];
            row.setAttribute("data-tt-parent-id", res[1]);
            if (typeof parent.children === "undefined"){
                _addClasse(parent.elem, 'selected');               
                parent.children = [];
            }
                
            parent.children.push(entry);
            if (rTestError.test(res[2])) {
                
                _addClasse(row, 'error-elem');
            }
        }
        collection[fullPath]= {
            elem: row
        };
    });
    $(table).treetable({
        expandable: true , 
        clickableNodeNames: true, 
        indent: 19, 
        indenterTemplate: '<span class="indenter" ><a href="#" title="Collapse">&nbsp;</a></span>'
    });
}




