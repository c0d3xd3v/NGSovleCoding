import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.15

Item {
    MouseArea{
        anchors.fill: parent
        onClicked: {
            if(!root_item.focus){
                root_item.focus = true
            }else if(root_item.focus)
            {
            }
        }
        ColumnLayout{
            anchors.fill: parent
            /*
            Item {
                Layout.fillWidth: true
                Layout.fillHeight: true
            }
            */
            GridLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                columns: 2
                Label{
                    text:"stop energy"
                }
                TextField{
                    id:stop_energy
                    Layout.fillWidth: true
                    placeholderText: "10"
                }
                Label{
                    text:"rel. edge length"
                }
                TextField{
                    id:rel_edge
                    Layout.fillWidth: true
                    placeholderText: "0.05"
                }
                Label{
                    text:"rel. epsilon"
                }
                TextField{
                    id: rel_eps
                    Layout.fillWidth: true
                    placeholderText: "0.001"
                }
            }
            Button {
                Layout.fillWidth: true
                text: "start meshing"
                onClicked: {
                    var se = parseFloat(stop_energy.text);
                    var rel = parseFloat(rel_edge.text);
                    var reps = parseFloat(rel_eps.text);
                    console.log(se)
                    mesh_controller.doMeshing(se, rel, reps);
                }
            }
            Item {
                Layout.fillWidth: true
                Layout.fillHeight: true
            }
        }
    }
}
