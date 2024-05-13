import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Item {
    id: item2
    x:10
    y:10
    width: 250
    height: 400
    visible: true
    clip: true

    property bool isFixed: false

    function updateFunctionList() {
        console.log("update")
        console.log(MainCtrl.getFunctionNames())
        comboBox.model = MainCtrl.getFunctionNames()
    }

    Component.onCompleted: function() {
        close()
        //MainCtrl.meshLoaded.connect(updateFunctionList)
    }

    function close(){
        if(!isFixed)
        {
            item2.anchors.top = undefined
            item2.anchors.bottom  = undefined
            height_range.to = roundButton.height
            width_range.to = roundButton.width
            widthAnim.start()
            heightAnim.start()
        }
    }

    //Component.onCompleted: close()

    SequentialAnimation on width {
        id: widthAnim
        running: false
        //PropertyAnimation { to: 36 }
        PropertyAnimation {id: width_range; to: 50 }
    }
    SequentialAnimation on height {
        id: heightAnim
        running: false
        //PropertyAnimation { to: 36 }
        PropertyAnimation {id: height_range; to: 50 }
        onFinished: function() {

            if(item2.height != roundButton.height)
            {
                item2.anchors.top = item2.parent.top
                item2.anchors.left  = item2.parent.left
                item2.anchors.bottom  = item2.parent.bottom
                item2.anchors.bottomMargin = 10
                item2.anchors.topMargin = 10
                item2.anchors.leftMargin = 10
            }
        }
    }

    Rectangle {
        color: "#82ffffff"
        anchors.fill: parent
        radius: 7
        border.color: "#fedcdcdc"
        border.width: 1

        MouseArea {
            id: mouseArea
            anchors.fill: parent
            anchors.topMargin: 0
            anchors.bottomMargin: 0
            anchors.leftMargin: 0
            anchors.rightMargin: 0
            focus: false
            hoverEnabled: true
            drag.threshold: 10
            onClicked: {
                if(width_range.to === 200){
                }else{
                    height_range.to = appWindow.height - 20
                    width_range.to = 200
                    widthAnim.start()
                    heightAnim.start()
                    //item2.anchors.bottom = appWindow.bottom
                }
                forceActiveFocus()
            }

            RowLayout {
                id: rowLayout
                height: 32
                anchors.left: parent.left
                anchors.right: parent.right
                anchors.top: parent.top
                anchors.topMargin: 0
                anchors.leftMargin: 0
                anchors.rightMargin: 0

                RoundButton {
                    id: roundButton1
                    opacity: 0.402
                    text: ""
                    icon.source: "../icons/settings.svg"
                    display: AbstractButton.IconOnly
                    onClicked: function() {
                        //item2.
                        height_range.to = appWindow.height - 20
                        width_range.to = 200
                        widthAnim.start()
                        heightAnim.start()
                    }
                }

                Item {
                    id: item3

                    Layout.fillHeight: true
                    Layout.fillWidth: true
                }


                RoundButton {
                    id: roundButton
                    opacity: 0.413
                    padding: 6
                    checked: false
                    checkable: true
                    highlighted: false
                    flat: false
                    display: AbstractButton.IconOnly
                    icon.source: "../icons/keep_off.svg"
                    onToggled: if(checked) {
                                   isFixed = true
                               } else if(!checked) {
                                   isFixed = false
                               }
                }



            }

            GridLayout {
                id: gridLayout
                anchors.left: parent.left
                anchors.right: parent.right
                anchors.top: rowLayout.bottom
                anchors.bottom: parent.bottom
                anchors.topMargin: 5
                columnSpacing: 5
                columns: 1

                RenderingControl {
                    id: columnLayout
                    height: 350
                }

                SolutionControl {
                    Layout.fillWidth: true
                    Layout.preferredHeight: 150
                }


                AnimationControls {
                    id: animationControls
                    Layout.fillWidth: true
                    height: 150
                }


                Item {
                    id: item1
                    Layout.minimumHeight: 150
                    Layout.fillHeight: true
                    Layout.fillWidth: true
                }


            }
        }

    }
}
