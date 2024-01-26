import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

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
        Component.onCompleted: {
             mesh_controller.onMeshingFinished.connect(function(){
                 var rep = mesh_controller.getRepresentation()
                 visCmb.currentIndex = visCmb.find("Result Mesh")
                 repCmb.currentIndex = repCmb.find(rep)
             })
        }

        ColumnLayout{
            anchors.fill: parent
            GridLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                columns: 2
                Label{
                    text:"visualization"
                }
                ComboBox{
                    id: visCmb
                    Layout.fillWidth: true
                    model: ["Source Mesh", "Result Mesh"]
                    onCurrentValueChanged: {
                        mesh_controller.setVisualization(currentValue)
                        var rep = mesh_controller.getRepresentation()
                        console.log("!!!!!!!!!!!!!!!! "  + rep)
                        var index = repCmb.find(rep)
                        repCmb.currentIndex = index
                    }
                }
                Label{
                    text:"representation"
                }
                ComboBox{
                    id: repCmb
                    Layout.fillWidth: true
                    model: ["Surface", "Wireframe", "Surface edges"]
                    onCurrentValueChanged: {
                        mesh_controller.setRepresentation(currentValue)
                    }
                }
            }
            Item {
                Layout.fillWidth: true
                Layout.fillHeight: true
            }
        }
    }
}
